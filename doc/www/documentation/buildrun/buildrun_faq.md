```{eval-rst}
:orphan:
```

(ref_buildrun_faq)=
# FAQ

## Configuring and Building

1. <a name="faq_confbuild1" href="#faq_confbuild1">**I run the configure script without any arguments and it fails. What should I do?**</a>

    First, we recommend checking whether there is a suitable [](ref_buildrun_configuration_wrappers) in the {{ '[`config`]({}/config)'.format(base_url) }} directory that you could use instead of running the configure script directly. If that is not the case, you need at least to specify the `LIBS` variable telling the configure script which libraries to link the executables to. The content of the list depends on the configure options you specify (see [Table 1](tab_icon_depgraph)), for example:

    ```sh
    ./configure --disable-mpi --disable-coupling LIBS='-lnetcdff -lnetcdf -llapack'
    ```

    If the libraries reside in nonstandard directories, you might also need to specify the `FCFLAGS`, `CPPFLAGS`, and `LDFLAGS` variables to tell the script which directories need to be searched for header and library files (see section [](ref_buildrun_configuration)) for more details).

2. <a name="faq_confbuild2" href="#faq_confbuild2">**How can I reproduce the configuration of the model used in a Buildbot test?**</a>

    Scripts run by Buildbot for configuration and building of the model reside in the `config/buildbot` directory. You can run them manually on the corresponding machine.

3. <a name="faq_confbuild3" href="#faq_confbuild3">**I get an error message from the configure script starting with _"configure: error: unable to find sources of..."_. What does this mean?**</a>

    Most probably, you forgot to initialize and/or update git submodules. You can do that by switching to the *source* root directory of ICON and running the following command:

    ```sh
    git submodule update --init
    ```

4. <a name="faq_confbuild4" href="#faq_confbuild4">**I have problems configuring/building ICON. What is the most efficient way to ask for help?**</a>

    Whoever you ask for help will appreciate receiving the log files. You can generate a tarball with the log files by running the following commands from the root build directory of ICON:

    ```sh
    make V=1 2>&1 | tee make.log
    tar --transform 's:^:build-report/:' -czf build-report.tar.gz $(find . -name 'config.log' -o -name 'CMakeCache.txt') make.log
    ```

    The result of the commands above will be file `build-report.tar.gz`, which should be attached to the very first email describing your problem. Please, do not forget to specify the **repository** and the **branch** that you experience the issue with, preferably in the form of a URL (e.g. {{ '[https://gitlab.dkrz.de/icon/icon-model/-/tree/icon-{}-public](https://gitlab.dkrz.de/icon/icon-model/-/tree/icon-{}-public)'.format(release, release) }}).


## Running

1. <a name="faq_run1" href="#faq_run1">**A simulation is not behaving in the way I expected, any advice to start troubleshooting?**</a>

    - Make sure that namelist names are correctly written in the runscript (ICON skips misspelled namelists and continues the run without exiting).


## Output

1. <a name="faq_output1" href="#faq_output1">**How can I find the available output variables of a configuration?**</a>

    Set `msg_level=15` (or larger) in the `&run_nml` namelist to print all available "ICON short names" of output variables for the current configuration in the logfile.
