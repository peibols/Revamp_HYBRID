#include "GpuHydroOpenCL.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <dlfcn.h>
#include <sstream>

#include "HydroProfile.h"

namespace {

using cl_int = int;
using cl_uint = unsigned int;
using cl_ulong = unsigned long;
using cl_bool = cl_uint;
using cl_bitfield = cl_ulong;
using cl_device_type = cl_bitfield;
using cl_mem_flags = cl_bitfield;
using cl_command_queue_properties = cl_bitfield;
using cl_context_properties = intptr_t;

using cl_platform_id = void*;
using cl_device_id = void*;
using cl_context = void*;
using cl_command_queue = void*;
using cl_program = void*;
using cl_kernel = void*;
using cl_mem = void*;
using cl_event = void*;

constexpr cl_int CL_SUCCESS = 0;
constexpr cl_bool CL_TRUE = 1;
constexpr cl_device_type CL_DEVICE_TYPE_GPU = 1u << 2;
constexpr cl_mem_flags CL_MEM_READ_ONLY = 1u << 2;
constexpr cl_mem_flags CL_MEM_WRITE_ONLY = 1u << 1;
constexpr cl_mem_flags CL_MEM_COPY_HOST_PTR = 1u << 5;
constexpr cl_uint CL_DEVICE_NAME = 0x102B;
constexpr cl_uint CL_PROGRAM_BUILD_LOG = 0x1183;

using PFN_clGetPlatformIDs = cl_int (*)(cl_uint, cl_platform_id*, cl_uint*);
using PFN_clGetDeviceIDs = cl_int (*)(cl_platform_id, cl_device_type, cl_uint, cl_device_id*, cl_uint*);
using PFN_clGetDeviceInfo = cl_int (*)(cl_device_id, cl_uint, size_t, void*, size_t*);
using PFN_clCreateContext = cl_context (*)(const cl_context_properties*, cl_uint, const cl_device_id*, void*, void*, cl_int*);
using PFN_clCreateCommandQueue = cl_command_queue (*)(cl_context, cl_device_id, cl_command_queue_properties, cl_int*);
using PFN_clCreateBuffer = cl_mem (*)(cl_context, cl_mem_flags, size_t, void*, cl_int*);
using PFN_clCreateProgramWithSource = cl_program (*)(cl_context, cl_uint, const char**, const size_t*, cl_int*);
using PFN_clBuildProgram = cl_int (*)(cl_program, cl_uint, const cl_device_id*, const char*, void*, void*);
using PFN_clGetProgramBuildInfo = cl_int (*)(cl_program, cl_device_id, cl_uint, size_t, void*, size_t*);
using PFN_clCreateKernel = cl_kernel (*)(cl_program, const char*, cl_int*);
using PFN_clSetKernelArg = cl_int (*)(cl_kernel, cl_uint, size_t, const void*);
using PFN_clEnqueueNDRangeKernel = cl_int (*)(cl_command_queue, cl_kernel, cl_uint, const size_t*, const size_t*, const size_t*, cl_uint, const cl_event*, cl_event*);
using PFN_clEnqueueReadBuffer = cl_int (*)(cl_command_queue, cl_mem, cl_bool, size_t, size_t, void*, cl_uint, const cl_event*, cl_event*);
using PFN_clFinish = cl_int (*)(cl_command_queue);
using PFN_clReleaseMemObject = cl_int (*)(cl_mem);
using PFN_clReleaseKernel = cl_int (*)(cl_kernel);
using PFN_clReleaseProgram = cl_int (*)(cl_program);
using PFN_clReleaseCommandQueue = cl_int (*)(cl_command_queue);
using PFN_clReleaseContext = cl_int (*)(cl_context);

struct OpenCLFns {
    PFN_clGetPlatformIDs clGetPlatformIDs = nullptr;
    PFN_clGetDeviceIDs clGetDeviceIDs = nullptr;
    PFN_clGetDeviceInfo clGetDeviceInfo = nullptr;
    PFN_clCreateContext clCreateContext = nullptr;
    PFN_clCreateCommandQueue clCreateCommandQueue = nullptr;
    PFN_clCreateBuffer clCreateBuffer = nullptr;
    PFN_clCreateProgramWithSource clCreateProgramWithSource = nullptr;
    PFN_clBuildProgram clBuildProgram = nullptr;
    PFN_clGetProgramBuildInfo clGetProgramBuildInfo = nullptr;
    PFN_clCreateKernel clCreateKernel = nullptr;
    PFN_clSetKernelArg clSetKernelArg = nullptr;
    PFN_clEnqueueNDRangeKernel clEnqueueNDRangeKernel = nullptr;
    PFN_clEnqueueReadBuffer clEnqueueReadBuffer = nullptr;
    PFN_clFinish clFinish = nullptr;
    PFN_clReleaseMemObject clReleaseMemObject = nullptr;
    PFN_clReleaseKernel clReleaseKernel = nullptr;
    PFN_clReleaseProgram clReleaseProgram = nullptr;
    PFN_clReleaseCommandQueue clReleaseCommandQueue = nullptr;
    PFN_clReleaseContext clReleaseContext = nullptr;
};

OpenCLFns &fns() {
    static OpenCLFns functions;
    return functions;
}

template <typename T>
bool loadSymbol(void *lib, const char *name, T &out) {
    out = reinterpret_cast<T>(dlsym(lib, name));
    return out != nullptr;
}

bool loadOpenCL(void *lib, std::string &error) {
    auto &cl = fns();
    bool ok = true;
    ok &= loadSymbol(lib, "clGetPlatformIDs", cl.clGetPlatformIDs);
    ok &= loadSymbol(lib, "clGetDeviceIDs", cl.clGetDeviceIDs);
    ok &= loadSymbol(lib, "clGetDeviceInfo", cl.clGetDeviceInfo);
    ok &= loadSymbol(lib, "clCreateContext", cl.clCreateContext);
    ok &= loadSymbol(lib, "clCreateCommandQueue", cl.clCreateCommandQueue);
    ok &= loadSymbol(lib, "clCreateBuffer", cl.clCreateBuffer);
    ok &= loadSymbol(lib, "clCreateProgramWithSource", cl.clCreateProgramWithSource);
    ok &= loadSymbol(lib, "clBuildProgram", cl.clBuildProgram);
    ok &= loadSymbol(lib, "clGetProgramBuildInfo", cl.clGetProgramBuildInfo);
    ok &= loadSymbol(lib, "clCreateKernel", cl.clCreateKernel);
    ok &= loadSymbol(lib, "clSetKernelArg", cl.clSetKernelArg);
    ok &= loadSymbol(lib, "clEnqueueNDRangeKernel", cl.clEnqueueNDRangeKernel);
    ok &= loadSymbol(lib, "clEnqueueReadBuffer", cl.clEnqueueReadBuffer);
    ok &= loadSymbol(lib, "clFinish", cl.clFinish);
    ok &= loadSymbol(lib, "clReleaseMemObject", cl.clReleaseMemObject);
    ok &= loadSymbol(lib, "clReleaseKernel", cl.clReleaseKernel);
    ok &= loadSymbol(lib, "clReleaseProgram", cl.clReleaseProgram);
    ok &= loadSymbol(lib, "clReleaseCommandQueue", cl.clReleaseCommandQueue);
    ok &= loadSymbol(lib, "clReleaseContext", cl.clReleaseContext);
    if (!ok) {
        error = "libOpenCL is present but required OpenCL 1.2 symbols are missing";
    }
    return ok;
}

std::string clError(const char *where, cl_int code) {
    std::ostringstream os;
    os << where << " failed with OpenCL error " << code;
    return os.str();
}

const char *kernelSource() {
    return R"CLC(
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

static int clamp_index(int idx, int max_count) {
    if (idx < 0) return 0;
    int max_idx = max_count > 1 ? max_count - 1 : 1;
    if (idx > max_idx - 1) return max_idx - 1;
    return idx;
}

static ulong hydro_index(int it, int iy, int ix, int ietamax, int ixmax) {
    return ((ulong)it * (ulong)ietamax + (ulong)iy) * (ulong)ixmax + (ulong)ix;
}

__kernel void hydro_interp(__global const double *tau_in,
                           __global const double *x_in,
                           __global const double *y_in,
                           __global double *temp_out,
                           __global double *vx_out,
                           __global double *vy_out,
                           __global const double *hydrot,
                           __global const double *hydrox,
                           __global const double *hydroy,
                           int n,
                           int itaumax,
                           int ietamax,
                           int ixmax,
                           double hydro_tau0,
                           double hydro_dtau,
                           double hydro_tau_max,
                           double hydro_xmax,
                           double hydro_dx,
                           double temp_scale) {
    int gid = (int)get_global_id(0);
    if (gid >= n) return;

    double tau = tau_in[gid];
    double x = x_in[gid];
    double y = y_in[gid];
    if (tau >= hydro_tau_max || tau < hydro_tau0) {
        temp_out[gid] = 0.0;
        vx_out[gid] = 0.0;
        vy_out[gid] = 0.0;
        return;
    }

    int it = clamp_index((int)floor((tau - hydro_tau0) / hydro_dtau), itaumax);
    double dt = (tau - (hydro_tau0 + (double)it * hydro_dtau)) / hydro_dtau;

    double xgrid = (hydro_xmax + x) / hydro_dx;
    int ix = clamp_index((int)floor(xgrid), ixmax);
    double dx = xgrid - (double)ix;

    double ygrid = (hydro_xmax + y) / hydro_dx;
    int iy = clamp_index((int)floor(ygrid), ietamax);
    double dy = ygrid - (double)iy;

    ulong base = hydro_index(it, iy, ix, ietamax, ixmax);
    ulong stride_y = (ulong)ixmax;
    ulong stride_t = (ulong)ietamax * (ulong)ixmax;

    double w000 = (1.0 - dt) * (1.0 - dx) * (1.0 - dy);
    double w100 = dt         * (1.0 - dx) * (1.0 - dy);
    double w010 = (1.0 - dt) * dx         * (1.0 - dy);
    double w001 = (1.0 - dt) * (1.0 - dx) * dy;
    double w110 = dt         * dx         * (1.0 - dy);
    double w011 = (1.0 - dt) * dx         * dy;
    double w101 = dt         * (1.0 - dx) * dy;
    double w111 = dt         * dx         * dy;

    ulong i000 = base;
    ulong i100 = base + stride_t;
    ulong i010 = base + 1;
    ulong i001 = base + stride_y;
    ulong i110 = base + stride_t + 1;
    ulong i011 = base + stride_y + 1;
    ulong i101 = base + stride_t + stride_y;
    ulong i111 = base + stride_t + stride_y + 1;

    temp_out[gid] = (hydrot[i000] * w000 + hydrot[i100] * w100 +
                     hydrot[i010] * w010 + hydrot[i001] * w001 +
                     hydrot[i110] * w110 + hydrot[i011] * w011 +
                     hydrot[i101] * w101 + hydrot[i111] * w111) * temp_scale;
    vx_out[gid] = hydrox[i000] * w000 + hydrox[i100] * w100 +
                  hydrox[i010] * w010 + hydrox[i001] * w001 +
                  hydrox[i110] * w110 + hydrox[i011] * w011 +
                  hydrox[i101] * w101 + hydrox[i111] * w111;
    vy_out[gid] = hydroy[i000] * w000 + hydroy[i100] * w100 +
                  hydroy[i010] * w010 + hydroy[i001] * w001 +
                  hydroy[i110] * w110 + hydroy[i011] * w011 +
                  hydroy[i101] * w101 + hydroy[i111] * w111;
}
)CLC";
}

}  // namespace

GpuHydroOpenCL::GpuHydroOpenCL() = default;

GpuHydroOpenCL::~GpuHydroOpenCL() {
    release();
}

void GpuHydroOpenCL::release() {
    auto &cl = fns();
    if (hydro_t_) cl.clReleaseMemObject(static_cast<cl_mem>(hydro_t_));
    if (hydro_vx_) cl.clReleaseMemObject(static_cast<cl_mem>(hydro_vx_));
    if (hydro_vy_) cl.clReleaseMemObject(static_cast<cl_mem>(hydro_vy_));
    if (kernel_) cl.clReleaseKernel(static_cast<cl_kernel>(kernel_));
    if (program_) cl.clReleaseProgram(static_cast<cl_program>(program_));
    if (queue_) cl.clReleaseCommandQueue(static_cast<cl_command_queue>(queue_));
    if (context_) cl.clReleaseContext(static_cast<cl_context>(context_));
    hydro_t_ = hydro_vx_ = hydro_vy_ = nullptr;
    kernel_ = program_ = queue_ = context_ = nullptr;
    platform_ = device_ = nullptr;
    if (opencl_lib_) {
        dlclose(opencl_lib_);
        opencl_lib_ = nullptr;
    }
}

bool GpuHydroOpenCL::initialize(const HydroProfile &hydro, std::string &error) {
    release();

    if (hydro.hydroT().empty() || hydro.hydroVx().empty() || hydro.hydroVy().empty()) {
        error = "hydro profile is empty";
        return false;
    }

    opencl_lib_ = dlopen("libOpenCL.so.1", RTLD_NOW);
    if (!opencl_lib_) opencl_lib_ = dlopen("libOpenCL.so", RTLD_NOW);
    if (!opencl_lib_) {
        error = "could not load libOpenCL";
        return false;
    }
    if (!loadOpenCL(opencl_lib_, error)) return false;

    auto &cl = fns();
    cl_uint platform_count = 0;
    cl_int err = cl.clGetPlatformIDs(0, nullptr, &platform_count);
    if (err != CL_SUCCESS || platform_count == 0) {
        error = clError("clGetPlatformIDs", err);
        return false;
    }

    std::vector<cl_platform_id> platforms(platform_count);
    err = cl.clGetPlatformIDs(platform_count, platforms.data(), nullptr);
    if (err != CL_SUCCESS) {
        error = clError("clGetPlatformIDs(list)", err);
        return false;
    }

    cl_device_id selected_device = nullptr;
    cl_platform_id selected_platform = nullptr;
    for (cl_platform_id platform : platforms) {
        cl_uint device_count = 0;
        err = cl.clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &device_count);
        if (err != CL_SUCCESS || device_count == 0) continue;
        std::vector<cl_device_id> devices(device_count);
        err = cl.clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, device_count, devices.data(), nullptr);
        if (err == CL_SUCCESS && !devices.empty()) {
            selected_platform = platform;
            selected_device = devices[0];
            break;
        }
    }
    if (!selected_device) {
        error = "no OpenCL GPU device found";
        return false;
    }

    platform_ = selected_platform;
    device_ = selected_device;

    size_t device_name_size = 0;
    cl.clGetDeviceInfo(selected_device, CL_DEVICE_NAME, 0, nullptr, &device_name_size);
    if (device_name_size > 0) {
        std::vector<char> name(device_name_size, '\0');
        cl.clGetDeviceInfo(selected_device, CL_DEVICE_NAME, name.size(), name.data(), nullptr);
        device_name_ = name.data();
    }

    cl_context context = cl.clCreateContext(nullptr, 1, &selected_device, nullptr, nullptr, &err);
    if (err != CL_SUCCESS || !context) {
        error = clError("clCreateContext", err);
        return false;
    }
    context_ = context;

    cl_command_queue queue = cl.clCreateCommandQueue(context, selected_device, 0, &err);
    if (err != CL_SUCCESS || !queue) {
        error = clError("clCreateCommandQueue", err);
        return false;
    }
    queue_ = queue;

    const char *source = kernelSource();
    size_t source_len = std::strlen(source);
    cl_program program = cl.clCreateProgramWithSource(context, 1, &source, &source_len, &err);
    if (err != CL_SUCCESS || !program) {
        error = clError("clCreateProgramWithSource", err);
        return false;
    }
    program_ = program;

    err = cl.clBuildProgram(program, 1, &selected_device, "", nullptr, nullptr);
    if (err != CL_SUCCESS) {
        size_t log_size = 0;
        cl.clGetProgramBuildInfo(program, selected_device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
        std::vector<char> log(log_size + 1, '\0');
        if (log_size > 0) {
            cl.clGetProgramBuildInfo(program, selected_device, CL_PROGRAM_BUILD_LOG, log_size, log.data(), nullptr);
        }
        std::ostringstream os;
        os << clError("clBuildProgram", err) << "\n" << log.data();
        error = os.str();
        return false;
    }

    cl_kernel kernel = cl.clCreateKernel(program, "hydro_interp", &err);
    if (err != CL_SUCCESS || !kernel) {
        error = clError("clCreateKernel", err);
        return false;
    }
    kernel_ = kernel;

    ixmax_ = hydro.ixmax();
    ietamax_ = hydro.ietamax();
    itaumax_ = hydro.itaumax();
    hydro_dx_ = hydro.hydroDx();
    hydro_xmax_ = hydro.hydroXmax();
    hydro_tau0_ = hydro.hydroTau0();
    hydro_dtau_ = hydro.hydroDtau();
    hydro_tau_max_ = hydro.hydroTauMax();
    temp_scaling_factor_ = hydro.tempScalingFactor();

    auto createReadOnly = [&](const std::vector<double> &data, const char *name) -> cl_mem {
        cl_int local_err = CL_SUCCESS;
        cl_mem buffer = cl.clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                          data.size() * sizeof(double),
                                          const_cast<double*>(data.data()), &local_err);
        if (local_err != CL_SUCCESS) {
            error = clError(name, local_err);
            return nullptr;
        }
        return buffer;
    };

    hydro_t_ = createReadOnly(hydro.hydroT(), "clCreateBuffer(hydroT)");
    if (!hydro_t_) return false;
    hydro_vx_ = createReadOnly(hydro.hydroVx(), "clCreateBuffer(hydroVx)");
    if (!hydro_vx_) return false;
    hydro_vy_ = createReadOnly(hydro.hydroVy(), "clCreateBuffer(hydroVy)");
    if (!hydro_vy_) return false;

    return true;
}

bool GpuHydroOpenCL::interpolate(const std::vector<double> &tau,
                                 const std::vector<double> &x,
                                 const std::vector<double> &y,
                                 std::vector<double> &temp,
                                 std::vector<double> &vx,
                                 std::vector<double> &vy,
                                 GpuHydroTiming *timing,
                                 std::string &error) {
    if (!context_ || !queue_ || !kernel_) {
        error = "OpenCL backend is not initialized";
        return false;
    }
    if (tau.size() != x.size() || tau.size() != y.size()) {
        error = "input vectors have inconsistent sizes";
        return false;
    }
    if (tau.empty()) {
        temp.clear();
        vx.clear();
        vy.clear();
        return true;
    }
    if (tau.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        error = "input vector is too large for kernel int indexing";
        return false;
    }

    auto &cl = fns();
    const auto start = std::chrono::steady_clock::now();
    const size_t bytes = tau.size() * sizeof(double);
    cl_int err = CL_SUCCESS;

    cl_mem tau_buf = cl.clCreateBuffer(static_cast<cl_context>(context_), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                       bytes, const_cast<double*>(tau.data()), &err);
    if (err != CL_SUCCESS) {
        error = clError("clCreateBuffer(tau)", err);
        return false;
    }
    cl_mem x_buf = cl.clCreateBuffer(static_cast<cl_context>(context_), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     bytes, const_cast<double*>(x.data()), &err);
    if (err != CL_SUCCESS) {
        cl.clReleaseMemObject(tau_buf);
        error = clError("clCreateBuffer(x)", err);
        return false;
    }
    cl_mem y_buf = cl.clCreateBuffer(static_cast<cl_context>(context_), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     bytes, const_cast<double*>(y.data()), &err);
    if (err != CL_SUCCESS) {
        cl.clReleaseMemObject(tau_buf);
        cl.clReleaseMemObject(x_buf);
        error = clError("clCreateBuffer(y)", err);
        return false;
    }

    cl_mem temp_buf = cl.clCreateBuffer(static_cast<cl_context>(context_), CL_MEM_WRITE_ONLY, bytes, nullptr, &err);
    if (err != CL_SUCCESS) {
        cl.clReleaseMemObject(tau_buf);
        cl.clReleaseMemObject(x_buf);
        cl.clReleaseMemObject(y_buf);
        error = clError("clCreateBuffer(temp)", err);
        return false;
    }
    cl_mem vx_buf = cl.clCreateBuffer(static_cast<cl_context>(context_), CL_MEM_WRITE_ONLY, bytes, nullptr, &err);
    if (err != CL_SUCCESS) {
        cl.clReleaseMemObject(tau_buf);
        cl.clReleaseMemObject(x_buf);
        cl.clReleaseMemObject(y_buf);
        cl.clReleaseMemObject(temp_buf);
        error = clError("clCreateBuffer(vx)", err);
        return false;
    }
    cl_mem vy_buf = cl.clCreateBuffer(static_cast<cl_context>(context_), CL_MEM_WRITE_ONLY, bytes, nullptr, &err);
    if (err != CL_SUCCESS) {
        cl.clReleaseMemObject(tau_buf);
        cl.clReleaseMemObject(x_buf);
        cl.clReleaseMemObject(y_buf);
        cl.clReleaseMemObject(temp_buf);
        cl.clReleaseMemObject(vx_buf);
        error = clError("clCreateBuffer(vy)", err);
        return false;
    }

    const int n = static_cast<int>(tau.size());
    cl_uint arg = 0;
    auto setArg = [&](size_t size, const void *value, const char *name) -> bool {
        cl_int local_err = cl.clSetKernelArg(static_cast<cl_kernel>(kernel_), arg++, size, value);
        if (local_err != CL_SUCCESS) {
            error = clError(name, local_err);
            return false;
        }
        return true;
    };

    bool ok = true;
    ok &= setArg(sizeof(cl_mem), &tau_buf, "clSetKernelArg(tau)");
    ok &= setArg(sizeof(cl_mem), &x_buf, "clSetKernelArg(x)");
    ok &= setArg(sizeof(cl_mem), &y_buf, "clSetKernelArg(y)");
    ok &= setArg(sizeof(cl_mem), &temp_buf, "clSetKernelArg(temp)");
    ok &= setArg(sizeof(cl_mem), &vx_buf, "clSetKernelArg(vx)");
    ok &= setArg(sizeof(cl_mem), &vy_buf, "clSetKernelArg(vy)");
    ok &= setArg(sizeof(cl_mem), &hydro_t_, "clSetKernelArg(hydroT)");
    ok &= setArg(sizeof(cl_mem), &hydro_vx_, "clSetKernelArg(hydroVx)");
    ok &= setArg(sizeof(cl_mem), &hydro_vy_, "clSetKernelArg(hydroVy)");
    ok &= setArg(sizeof(int), &n, "clSetKernelArg(n)");
    ok &= setArg(sizeof(int), &itaumax_, "clSetKernelArg(itaumax)");
    ok &= setArg(sizeof(int), &ietamax_, "clSetKernelArg(ietamax)");
    ok &= setArg(sizeof(int), &ixmax_, "clSetKernelArg(ixmax)");
    ok &= setArg(sizeof(double), &hydro_tau0_, "clSetKernelArg(tau0)");
    ok &= setArg(sizeof(double), &hydro_dtau_, "clSetKernelArg(dtau)");
    ok &= setArg(sizeof(double), &hydro_tau_max_, "clSetKernelArg(tauMax)");
    ok &= setArg(sizeof(double), &hydro_xmax_, "clSetKernelArg(xmax)");
    ok &= setArg(sizeof(double), &hydro_dx_, "clSetKernelArg(dx)");
    ok &= setArg(sizeof(double), &temp_scaling_factor_, "clSetKernelArg(tempScale)");
    if (!ok) {
        cl.clReleaseMemObject(tau_buf);
        cl.clReleaseMemObject(x_buf);
        cl.clReleaseMemObject(y_buf);
        cl.clReleaseMemObject(temp_buf);
        cl.clReleaseMemObject(vx_buf);
        cl.clReleaseMemObject(vy_buf);
        return false;
    }

    const size_t local = 128;
    const size_t global = ((tau.size() + local - 1) / local) * local;
    err = cl.clEnqueueNDRangeKernel(static_cast<cl_command_queue>(queue_), static_cast<cl_kernel>(kernel_),
                                    1, nullptr, &global, &local, 0, nullptr, nullptr);
    if (err != CL_SUCCESS) {
        cl.clReleaseMemObject(tau_buf);
        cl.clReleaseMemObject(x_buf);
        cl.clReleaseMemObject(y_buf);
        cl.clReleaseMemObject(temp_buf);
        cl.clReleaseMemObject(vx_buf);
        cl.clReleaseMemObject(vy_buf);
        error = clError("clEnqueueNDRangeKernel", err);
        return false;
    }

    temp.resize(tau.size());
    vx.resize(tau.size());
    vy.resize(tau.size());
    err = cl.clEnqueueReadBuffer(static_cast<cl_command_queue>(queue_), temp_buf, CL_TRUE, 0, bytes,
                                 temp.data(), 0, nullptr, nullptr);
    if (err == CL_SUCCESS) {
        err = cl.clEnqueueReadBuffer(static_cast<cl_command_queue>(queue_), vx_buf, CL_TRUE, 0, bytes,
                                     vx.data(), 0, nullptr, nullptr);
    }
    if (err == CL_SUCCESS) {
        err = cl.clEnqueueReadBuffer(static_cast<cl_command_queue>(queue_), vy_buf, CL_TRUE, 0, bytes,
                                     vy.data(), 0, nullptr, nullptr);
    }
    if (err == CL_SUCCESS) {
        err = cl.clFinish(static_cast<cl_command_queue>(queue_));
    }

    cl.clReleaseMemObject(tau_buf);
    cl.clReleaseMemObject(x_buf);
    cl.clReleaseMemObject(y_buf);
    cl.clReleaseMemObject(temp_buf);
    cl.clReleaseMemObject(vx_buf);
    cl.clReleaseMemObject(vy_buf);

    if (err != CL_SUCCESS) {
        error = clError("readback/finish", err);
        return false;
    }

    if (timing) {
        const auto stop = std::chrono::steady_clock::now();
        timing->milliseconds = std::chrono::duration<double, std::milli>(stop - start).count();
        timing->device_name = device_name_;
    }
    return true;
}
