(ref_development_guidelines)=
# Development Guidelines

Welcome to the development guidelines of the ICON model!
For a smooth experience for all involved parties, please have a look at the **[Development Workflow](ref_dev_overview)** described below.
This contains also guidelines for [contributors](ref_contribute_guidelines), [reviewers](ref_review_guidelines) and [gatekeepers](ref_gatekeeping_guidelines).

A code contribution should always be accompanied by a corresponding documentation for [docs.icon-model.org](https://docs.icon-model.org).
The **[Documentation Guidelines](ref_docs_guidelines)** assist you in the generation of documentation for ICON.
Besides a description of the [publication workflow](ref_docs_guidelines_workflow), it provides you with a [starting point](ref_docs_guidelines_howtostart) including best practices and many [examples](ref_docs_guidelines_examples) for specific layout options you might want to use.

```{button-ref} ref_docs_guidelines
:color: primary
:outline:
To the Documentation Guidelines
```


(ref_dev_overview)=
## Development Workflow
{{author}}Ralf Mueller{{endauthor}}

This document explains the internal ICON development workflow: How do _contributors_
get code changes into the main development branches and the public releases.
The main tool in this process is the reviewing of the respective merge requestes.
Although the final merge will be performed by a comparably small group of
_gatekeepers_, the main resposibility lies with the _reviewers_.

This diagram illustrates the development workflow and specifies the roles of the contributor,
the reviewer and the gatekeeper.

```{image} wf.svg
:align: right
:height: 840
:width: 600
:alt: ICON development workflow
:name: fig_icon_workflow
```

Please note that the main ICON development is currently done in two major
groups:

- DWD/KIT/MCH: numerical weather and seasonal climate predictions as well as climate projections performed by DWD and other weather services
- MPI-MET: climates physics, dynamics, variability and related research

Both groups have slightly different development workflows, which will be referenced in
the following text and documents wherever necessary.

Depending on your role, check the following documents:
- [Contributor Guidelines](ref_contribute_guidelines)
- [Reviewer Guidelines](ref_review_guidelines)
- [Gatekeeper Guidelines](ref_gatekeeping_guidelines)

### Prerequisites for successful merges
(_for all contributors_)

- **Only one feature per merge request is allowed!** <br> Large
  features/changes exceeding a few hundreds of code lines should preferably be
  split into several steps. Existing
  {{ '[templates]({}/.gitlab/merge_request_templates)'.format(base_url) }}
  must be used.
- **Short, but precise description of the changeset and scope of the merge
  request** <br> For larger changes, more details should be provided in an
  accompanying GitLab issue
- **The feature branch must be rebased on the current state of the target branch**
- **Successful testing**
  - MPI-MET: All shown CI pipelines for the merge request (MR) have to be executed sucessfully
  - DWD:
    1. Buildbot Tests have to be green <br> In case of an outage of
       certain builders, the Gatekeepers can decide to ignore these builders.
    1. Mandatory tests for the last commit: `runBB` + manual run of CSCS-CI
- **The merge request must have undergone a full review** <br> Selecting reviewers is a contributor's task. Please consider
  the skills needed for a proper review of all introduced changes. Multiple
  reviewers might be needed. If in doubt ask gatekeepers for suggestions.
  - Reviewers must be added to the merge request in that role
  - Merge request is considered to be under review as long as a reviewer is
  assigned, but has not approved, yet<br>For DWD: The contributor is recommended to add the GitLab-label **`inReview`** until the reviewers approve the merge request.
  - **After successful Review: Reviewer uses Approve button!**
- **Development finished:**  contributors mark MR as **ready** (Draft = No) and unassign. This will schedule the merge request for getting a gatekeeper assigned.

When a merge request does not meet the above criteria, it will be reset to the `Draft`
status by the gatekeeper group. For a more detailed guide on how to work on your feature, please
read the [Contributor guidelines](ref_contribute_guidelines) and check the [Reviewer Guidelines](ref_review_guidelines).

```{button-ref} ref_contribute_guidelines
:color: primary
:outline:
To the Contribution Guidelines
```

```{button-ref} ref_review_guidelines
:color: primary
:outline:
To the Reviewing Guidelines
```

### Gatekeeping

Gatekeepers manage the merge process and make the final decision about whether
the code enters the master/main branch according to the
[Gatekeeper Guidelines](ref_gatekeeping_guidelines). They can request additional reviews or
further changes if needed.

- The group of gatekeepers will decide on the **assignee**
  - DWD: in a regular meeting (usually Mondays)
  - MPI-MET: in an anarchic collective decision process (usually within 3-5 days)
- **Resolving issues, monitoring tests, updating the reference data and rebasing stays in the responsibility of the contributor**.
- The assigned gatekeeper decides whether retriggering tests after the final
  rebase is necessary.

In the event that the gatekeepers are unable to reach a decision on whether to accept
or reject a merge request (even after following further review steps), the decision shall
be made by the C5. For more information, please refer to ICON's [governance model](https://www.icon-model.org/about-us/icon-governance)

```{button-ref} ref_gatekeeping_guidelines
:color: primary
:outline:
To the Gatekeeping Guidelines
```
