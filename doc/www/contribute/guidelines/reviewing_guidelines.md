```{eval-rst}
:orphan:
```

(ref_review_guidelines)=
# Reviewer Guidelines
{{author}}Ralf Mueller{{endauthor}}


The main task of the reviewer is to ensure the four-eyes principle. Reviewers
should be aware of the fact that no additional in-depth review will follow their
assessment. The final gatekeeper action only covers formal aspects.

Reviewers should be aware of the their own level of expertise. Changes that
go beyond shall be reviewed by another person. If in doubt comment on the
MR regarding the changes that have been reviewed by you so that the contributor is
aware that there is need for an additional reviewer.

The responsibilities of the reviewer are **Correctness**, **Performance**,
**Documentation** and **Appearance**.

## Correctness

- All exsting tests must work
- New code must be covered by tests
- Reviewers verify that the implementation is correct and the results are scientifically sound
- New code must be fully and correctly parallelized with MPI and OpenMP
- New code for operational NWP:
  - must vectorize
  - must work in combination with nesting (at least technically)
- GPU ports must follow OpenACC guidelines[^1]
  - If the changes should not run on GPU, they must invoke `CALL finish` as a
    safeguard
- Code design:
  - Is the new code well designed, e.g. encapsulation, proper memory layout, reuse of existing code infrastructure?
  - Unnecessary code duplication shall be avoided and the code shall be future proof
  - Anything that changes established defaults must be very well justified and
    agreed upon
  - Avoid unwanted side-effects of the added changes

## Performance

- Timers and memory consumption should be neutral for operationally used code
  parts (except for well-justified, agreed changes). The contributor has to
  provide proof in doubt (e.g. ICON's timer output, ...)
- Performance considerations during the setup phase can be handled less strict,
  particularly with respect to OpenMP parallelization.
- Constant operations like namelist cross-checks have to be done during the
  setup phase, not within the time loop.

## Documentation

The merge request (MR) description is the most important documentation outside the
repository. Its title and decription has to be in-line with the changes added
and as clearly as possible.

The connection between the underlying physics and the source code must be well
documented.  More specifically, this includes the following:

### In-code
- Variable declarations come with a description of the variable and the unit (SI!)
- No temporary change of the variable content (e.g. changing from a mass to a volume)
- Relevant literature should be referenced in a comment
### External
- Public documentation under `doc/www` must be adapted for new features and fixes
- Changes to namelist switches have to be properly documented in the namelist
  tables and the merge request description
- The Merge Request contains an overview on the changes and the
  motivation/reason for the changes
- If there is a change in a git submodule, these changes that come with the
  submodule update should be summarized and, if possible, the corresponding
  merge request(s) in the submodule repository referenced
### License/Maintenance
- Open Source: Code has to follow ICON's open source regulations (file header,
  third-party code has to come with a proper open source license, no author
  information in source files)
- Is it clear who will be responsible for the newly introduced code parts and
  the test scripts in the future? This should ideally take into account
  changesets created by non-permanent staff members.

## Appearance

New code parts should blend into the existing code, if possible even improve the
code readability.  As a rule of thumb, the coding style of new code should be
similar to other code in the corresponding module. See [^1] for basic guidelines.

Extensive code formatting should be done in dedicated merge requests.

Note: This will be done automatically by code formaters in the future

[^1]: _Further information: [https://gitlab.dkrz.de/icon/wiki/-/wikis/GPU-development/ICON-OpenAcc-style-and-implementation-guide](https://gitlab.dkrz.de/icon/wiki/-/wikis/GPU-development/ICON-OpenAcc-style-and-implementation-guide)_
