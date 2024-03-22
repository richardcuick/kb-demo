jupyter-book build %1\.
md dist\%1
xcopy %1\_build\html\*.* dist\%1 /E
