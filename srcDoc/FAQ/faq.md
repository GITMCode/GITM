# Common Problems

## The code won't compile!!

Sorry about that. It is hard to test on every system type. 

There are some steps you can try before submitting a bug report or contacting the developers.

1. First, try to un-install everything:
```
make distclean
./Config.pl -uninstall
```
2. Then, check to see if any files are changed:
```
git status
```
3. If any of the Fortran or Make-files have changed, you may want to `git restore` the changes.
4. Try to re-configure and then re-compile

