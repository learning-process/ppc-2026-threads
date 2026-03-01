# Run clang-tidy on sabutay_sparse_complex_ccs_mult sources.
# Google Test: C:\deps\gtest\include (copied from 3rdparty). Uses 3rdparty for full run.
# Run from project root: .\tasks\sabutay_sparse_complex_ccs_mult\run_clang_tidy.ps1

$ProjectRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$TaskDir = Join-Path $ProjectRoot "tasks\sabutay_sparse_complex_ccs_mult"

$IncludeArgs = @(
    "-I$ProjectRoot",
    "-I$ProjectRoot\tasks",
    "-I$ProjectRoot\modules",
    "-I$ProjectRoot\3rdparty",
    "-I$ProjectRoot\3rdparty\googletest\googletest\include",
    "-I$ProjectRoot\3rdparty\json\include",
    "-I$ProjectRoot\3rdparty\libenvpp\include",
    "-I$ProjectRoot\3rdparty\libenvpp\external\fmt\include",
    "-I$ProjectRoot\3rdparty\onetbb\include",
    "-std=c++20",
    "-DPPC_SETTINGS_sabutay_sparse_complex_ccs_mult=`"$TaskDir\settings.json`""
)

Set-Location $ProjectRoot
Get-ChildItem -Recurse -Filter *.cpp -Path $TaskDir | ForEach-Object {
    Write-Host "Checking $($_.FullName)..."
    & clang-tidy $_.FullName --checks="-modernize-use-scoped-lock" -- @IncludeArgs
}
