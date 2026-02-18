# C++ Conversion Build

Preferred (CMake), from `cpp conversion`:

```powershell
cmake -S . -B build
cmake --build build --config Release
```

Executable output:

- `build/ogcode_cpp.exe` (single-config generators)
- or `build/Release/ogcode_cpp.exe` (multi-config generators)

Fallback (direct `g++`), from project root:

```powershell
$src = Get-ChildItem "cpp conversion/src/*.cpp" | ForEach-Object { $_.FullName }
g++ -std=c++17 -I"cpp conversion/include" $src -o "cpp conversion/ogcode_cpp.exe"
```

If your toolchain supports OpenMP, add `-fopenmp`.

Run from the original project root (or copy required input files into the run directory) so the executable can find:

- `A_InputParameters.txt`
- and any other runtime text files used by the model.
