#pragma once

#include <cstddef>
#include <string>
#include <vector>

class HydroProfile;

struct GpuHydroTiming {
    double milliseconds = 0.0;
    std::string device_name;
};

class GpuHydroOpenCL {
public:
    GpuHydroOpenCL();
    ~GpuHydroOpenCL();

    GpuHydroOpenCL(const GpuHydroOpenCL&) = delete;
    GpuHydroOpenCL& operator=(const GpuHydroOpenCL&) = delete;

    bool initialize(const HydroProfile &hydro, std::string &error);
    bool interpolate(const std::vector<double> &tau,
                     const std::vector<double> &x,
                     const std::vector<double> &y,
                     std::vector<double> &temp,
                     std::vector<double> &vx,
                     std::vector<double> &vy,
                     GpuHydroTiming *timing,
                     std::string &error);

private:
    void release();

    void *opencl_lib_ = nullptr;
    void *platform_ = nullptr;
    void *device_ = nullptr;
    void *context_ = nullptr;
    void *queue_ = nullptr;
    void *program_ = nullptr;
    void *kernel_ = nullptr;
    void *hydro_t_ = nullptr;
    void *hydro_vx_ = nullptr;
    void *hydro_vy_ = nullptr;

    int ixmax_ = 0;
    int ietamax_ = 0;
    int itaumax_ = 0;
    double hydro_dx_ = 0.0;
    double hydro_xmax_ = 0.0;
    double hydro_tau0_ = 0.0;
    double hydro_dtau_ = 0.0;
    double hydro_tau_max_ = 0.0;
    double temp_scaling_factor_ = 1.0;
    std::string device_name_;
};
