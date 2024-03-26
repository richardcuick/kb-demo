set PYDEVD_DISABLE_FILE_VALIDATION=1
rmdir /s /q %1\_build
rmdir /s /q dist\%1
jupyter-book build %1\.
md dist\%1
xcopy %1\_build\html\*.* dist\%1 /E /Y
