## Changes in the new version

The code is essentially the same, just some simple cleanup and added a new makefile (`makefile.vc`) for compining in Windows machines:

```
nmake /f makefile.vc clean
nmake /f makefile.vc
```

Please note that the `-mode 4` option is no longer available - it was for research only and no longer needed.

