# General information

This directory contains a collection of common building blocks (*templates*) for ICON GitLab CI jobs. They are implemented as [hidden jobs](https://docs.gitlab.com/ci/jobs/#hide-a-job), which can be [extended](https://docs.gitlab.com/ci/yaml/#extends) and [referenced](https://docs.gitlab.com/ci/yaml/yaml_optimization/#reference-tags) by the normal ones. The templates in this directory are fully functional out-of-the-box (i.e. a normal job that implements a [`script`](https://docs.gitlab.com/ci/yaml/#script) and `extends` one of the templates is supposed to be picked up and run by a GitLab runner).

# Requirements and conventions

1. Names of all jobs defined in `.gitlab/ci/common/<filename>.yml` (i.e. in this directory) must be `.common:<filename>` or start with `.common:<filename>:`. This makes it easy to find the file with the definition by the job name.
2. The jobs are not allowed to specify any [job keywords](https://docs.gitlab.com/ci/yaml/#job-keywords) that can affect the pipeline of the extending jobs (e.g. [`allow_failure`](https://docs.gitlab.com/ci/yaml/#allow_failure), [`interruptible`](https://docs.gitlab.com/ci/yaml/#interruptible), [`rules`](https://docs.gitlab.com/ci/yaml/#rules), [`when`](https://docs.gitlab.com/ci/yaml/#when), etc.). Exceptions to this rule are possible but must be discussed during the review.
3. It is recommended to extend the [`.utils:log`](/.gitlab/ci/utils/log.yml) template and implement the [`before_script`](https://docs.gitlab.com/ci/yaml/#before_script) of the job as a single multi-line string as follows:
    ```yaml
    include:
      - local: .gitlab/ci/utils/log.yml

    .common:whatever:job:
      extends: .utils:log
      before_script:
        - !reference [.utils:log, before_script]
        - |
          # !reference[.common:whatever:job, before_script]
          icon_ci_log_section_start 'Doing whatever'
          (
            whatever_var='whatever value'
            set -ux
            whatever_cmd "${whatever_var}"
          ) 2>&1
          icon_ci_log_section_end
    ```
    This way, the log of the job that extends the template is cleaner and more informative for advanced users. The benefits of implementing the script as single subshell are listed in [`.gitlab/ci/utils/README.md`](/.gitlab/ci/utils/README.md).
4. Names of all custom [`variables`](https://docs.gitlab.com/ci/yaml/#job-variables) defined by a job must start with `ICON_CI_COMMON_<FILENAME>_`. The `ICON_` prefix tells that it is a custom variable, which is not supposed to affect any tools run by the job. The `CI_` infix tells that the variable is CI-specific and is supposed to affect neither ICON scripts nor ICON itself. The `COMMON_<FILENAME>_` infix makes it easy to find the file with the definition of the job by the variable name.
