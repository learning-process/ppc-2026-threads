{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  nativeBuildInputs = [
    pkgs.gcc
    pkgs.gdb
    pkgs.cmake
    pkgs.pkg-config
    pkgs.ninja
    pkgs.openmpi
  ];

  buildInputs = [
    pkgs.openmpi
  ];

  shellHook = ''
    export CC=${pkgs.gcc}/bin/gcc
    export CXX=${pkgs.gcc}/bin/g++
    export OMPI_CC=${pkgs.gcc}/bin/gcc
    export OMPI_CXX=${pkgs.gcc}/bin/g++

    # Пути к локальным библиотекам (после сборки)
    export LD_LIBRARY_PATH="$PWD/build/ppc_googletest/install/lib64:$PWD/build/ppc_libenvpp/install/lib64:$PWD/build/ppc_onetbb/install/lib64:$LD_LIBRARY_PATH"
    export LIBRARY_PATH="$LD_LIBRARY_PATH"
    export CPLUS_INCLUDE_PATH="$PWD/build/ppc_googletest/install/include:$PWD/build/ppc_libenvpp/install/include:$PWD/build/ppc_onetbb/install/include:$CPLUS_INCLUDE_PATH"
    echo "Environment setup for PPC (using local 3rdparty libraries)"
  '';
}
