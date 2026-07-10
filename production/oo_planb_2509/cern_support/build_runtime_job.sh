#!/usr/bin/env bash
set -Eeuo pipefail
EOS_BASE="${EOS_BASE:?EOS_BASE is required}"
WORK="$PWD/work"
INITIAL_DIR="$PWD"
fetch() { xrdcp -f "root://eosuser.cern.ch/${EOS_BASE}/$1" "$2"; }
put() { xrdcp -f "$1" "root://eosuser.cern.ch/${EOS_BASE}/$2"; }
finalize() {
  rc=$?; set +e
  echo "status=$([[ $rc -eq 0 ]] && echo success || echo failed)" > "$INITIAL_DIR/build_status.txt"
  echo "exit_code=$rc" >> "$INITIAL_DIR/build_status.txt"
  echo "date=$(date -Is)" >> "$INITIAL_DIR/build_status.txt"
  put "$INITIAL_DIR/build_status.txt" status/build_status.txt >/dev/null 2>&1 || true
  exit $rc
}
trap finalize EXIT
mkdir -p "$WORK"
cd "$WORK"
fetch payloads/pythia8315_alma9_install.tar.gz pythia8315_alma9_install.tar.gz
fetch payloads/mmli_source.tar.gz mmli_source.tar.gz
fetch payloads/runtime_payload.tar.gz runtime_payload.tar.gz
tar -xzf pythia8315_alma9_install.tar.gz
tar -xzf runtime_payload.tar.gz
mkdir -p source/mmli bin
tar -xzf mmli_source.tar.gz -C source/mmli
PYTHIA_ROOT="$WORK/pythia8315_alma9_install"
export PYTHIA_INCLUDE="$PYTHIA_ROOT/include"
export PYTHIA_LIB="$PYTHIA_ROOT/lib"
export PYTHIA8="$PYTHIA_ROOT"
export PYTHIA8DATA="$PYTHIA_ROOT/share/Pythia8/xmldoc"
export LD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${LD_LIBRARY_PATH:-}"
ROOT_LIBDIR="${ROOT_LIBDIR:-}"
if [[ -z "$ROOT_LIBDIR" ]] && command -v root-config >/dev/null 2>&1; then
  ROOT_LIBDIR="$(root-config --libdir)"
fi
if [[ -z "$ROOT_LIBDIR" && -d /usr/lib64/root ]]; then
  ROOT_LIBDIR="/usr/lib64/root"
fi
if [[ -n "$ROOT_LIBDIR" ]]; then
  test -d "$ROOT_LIBDIR"
  echo "$ROOT_LIBDIR" > "$WORK/runtime/root_libdir.txt"
  export LD_LIBRARY_PATH="$ROOT_LIBDIR:$LD_LIBRARY_PATH"
fi
LHAPDF_VIEW="${LHAPDF_CVMFS_VIEW:-}"
if [[ -z "$LHAPDF_VIEW" && -f "$WORK/runtime/lhapdf_cvmfs_view.txt" ]]; then
  LHAPDF_VIEW="$(cat "$WORK/runtime/lhapdf_cvmfs_view.txt")"
fi
if [[ -n "$LHAPDF_VIEW" ]]; then
  test -d "$LHAPDF_VIEW/include/LHAPDF"
  test -f "$LHAPDF_VIEW/lib/libLHAPDF.so"
  export LD_LIBRARY_PATH="$LHAPDF_VIEW/lib:$LD_LIBRARY_PATH"
fi
if [[ -d "$WORK/runtime/lhapdf_data" ]]; then
  export LHAPDF_DATA_PATH="$WORK/runtime/lhapdf_data${LHAPDF_DATA_PATH:+:$LHAPDF_DATA_PATH}"
fi
if [[ -n "$LHAPDF_VIEW" && -d "$LHAPDF_VIEW/share/LHAPDF" ]]; then
  export LHAPDF_DATA_PATH="${LHAPDF_DATA_PATH:+$LHAPDF_DATA_PATH:}$LHAPDF_VIEW/share/LHAPDF"
fi
if [[ -n "$LHAPDF_VIEW" ]]; then
  mkdir -p "$WORK/runtime/lhapdf_plugin"
  cat > "$WORK/LHAPDF6Plugin.cc" <<'EOF_PLUGIN'
#include "Pythia8Plugins/LHAPDF6.h"
EOF_PLUGIN
  g++ -std=c++17 -O2 -fPIC -shared \
    -I"$PYTHIA_ROOT/include" -I"$LHAPDF_VIEW/include" \
    "$WORK/LHAPDF6Plugin.cc" \
    -L"$PYTHIA_ROOT/lib" -L"$LHAPDF_VIEW/lib" \
    -Wl,-rpath,"$LHAPDF_VIEW/lib" -Wl,-rpath,"$PYTHIA_ROOT/lib" \
    -lLHAPDF -lpythia8 \
    -o "$WORK/runtime/lhapdf_plugin/libpythia8lhapdf6.so" \
    > "$WORK/lhapdf6_plugin_build.log" 2>&1
fi
(
  cd source/mmli
  ./compiler.sh main > "$WORK/mmli_build.log" 2>&1
  cp main "$WORK/bin/main"
)
ldd bin/main > mmli_main_ldd.txt
if [[ -f "$WORK/runtime/lhapdf_plugin/libpythia8lhapdf6.so" ]]; then
  ldd "$WORK/runtime/lhapdf_plugin/libpythia8lhapdf6.so" > lhapdf6_plugin_ldd.txt || true
fi
bundle_root_libs() {
  local bundle_dir="$WORK/runtime/root_lib"
  local ldd_file soname libpath resolved realbase
  : > "$WORK/root_bundle_manifest.txt"
  : > "$WORK/root_bundle_size.txt"
  if [[ -z "${ROOT_LIBDIR:-}" || ! -d "$ROOT_LIBDIR" ]]; then
    echo "root_libdir=none" >> "$WORK/root_bundle_manifest.txt"
    return 0
  fi
  mkdir -p "$bundle_dir"
  echo "root_libdir=$ROOT_LIBDIR" >> "$WORK/root_bundle_manifest.txt"
  for ldd_file in "$WORK/mmli_main_ldd.txt" "$WORK/lhapdf6_plugin_ldd.txt"; do
    [[ -f "$ldd_file" ]] || continue
    while read -r soname libpath; do
      [[ -n "$soname" && -f "$libpath" ]] || continue
      case "$libpath" in
        "$ROOT_LIBDIR"/*)
          resolved="$(readlink -f "$libpath")"
          realbase="$(basename "$resolved")"
          if [[ ! -f "$bundle_dir/$realbase" ]]; then
            install -m 755 "$resolved" "$bundle_dir/$realbase"
          fi
          if [[ "$soname" != "$realbase" && ! -e "$bundle_dir/$soname" ]]; then
            ln -s "$realbase" "$bundle_dir/$soname"
          fi
          printf '%s\t%s\t%s\n' "$soname" "$libpath" "$resolved" >> "$WORK/root_bundle_manifest.txt"
          ;;
      esac
    done < <(awk '/=>/ {print $1, $(NF-1)} /^[[:space:]]*\/.*\.so/ {print $1, $1}' "$ldd_file")
  done
  find "$bundle_dir" -maxdepth 1 -name '*.so*' -printf 'bundled\t%f\n' | sort >> "$WORK/root_bundle_manifest.txt" || true
  du -sh "$bundle_dir" > "$WORK/root_bundle_size.txt" || true
}
bundle_root_libs
tar -czf mmli_runtime_alma9.tar.gz bin runtime mmli_build.log mmli_main_ldd.txt lhapdf6_plugin_build.log lhapdf6_plugin_ldd.txt root_bundle_manifest.txt root_bundle_size.txt
put mmli_runtime_alma9.tar.gz payloads/mmli_runtime_alma9.tar.gz
