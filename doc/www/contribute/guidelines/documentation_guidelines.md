```{eval-rst}
:orphan:
```

(ref_docs_guidelines)=
# Documentation Guidelines
{{author}}D. Rieger{{endauthor}}

We want to encourage all ICON developers in contributing to the [ICON documentation (docs.icon-model.org)](https://docs.icon-model.org).
Your contributions to improve the coverage and readability of this website are necessary and always welcome!

The ICON documentation is written in an extended form of Markdown, which should be a familiar format for users of Gitlab, Github and thus ease the creation of documentation.
Further information on the software used to generate the HTML website from Markdown is available on the websites listed in the [useful links section](ref_docs_guidelines_links).

The following section provide you with some information on the [publication workflow](ref_docs_guidelines_workflow), a [starting point](ref_docs_guidelines_howtostart), the [terms of reference](ref_docs_guidelines_tor) and [several examples](ref_docs_guidelines_examples).

(ref_docs_guidelines_workflow)=
## Publication Workflow

The documentation is located in the subfolder {{ '[`doc/www`]({}/doc/www)'.format(base_url) }}.
Changes in this subfolder will trigger the CI/CD pipeline and generate a HTML preview.
A **View app** button will appear in the merge request overview website (right below the pipelines) which will lead you to the HTML preview.

When a new ICON release is published, the [docs.icon-model.org Website](https://docs.icon-model.org) is deployed correspondingly.
A preview of the content from the current `icon:main` is deployed under [https://icon.gitlab-pages.dkrz.de/icon](https://icon.gitlab-pages.dkrz.de/icon).
Both are publicly accessible.

(ref_docs_guidelines_howtostart)=
## How to Start

Since the documentation is part of the ICON repository, documentation updates take place in merge requests, following the same workflow as code developments.
Ideally, code developments should always be accompanied by documentation updates.
Merge request can be used to experiment with modifications in the {{ '[`doc/www`]({}/doc/www)'.format(base_url) }} subfolder.
There exist already many solutions and examples for specific formats, admonition, tables or figures, some are outlined in detail [below](ref_docs_guidelines_examples).
Probably you will find an example/solution for your problem somewhere in the code or in the [documentation of the software used to generate the website](ref_docs_guidelines_links).

:::{admonition} Best Practices
:class: admonition-icontheme
* Please settle all data protection issues in advance.
* Please keep the number of warnings at zero. The pipeline will fail otherwise.
* As there is no automatic mail obfuscation in Sphinx, you might want to obscure mail adresses manually to avoid Spam
* Please locate static content like images and pdfs used in the documentation in the same directory as the file where it is referenced to ease (relative) referencing, moving and tidying up content.
* New subsites should not be added to the `toctree` in {{ '[`index.md`]({}/doc/www/index.md)'.format(base_url) }} since this information is used to generate the links in the header of the website.
:::

(ref_docs_guidelines_examples)=
## Examples

### Tables

A table with a caption can be created like this:

:::{table} Table 1: Example table
:width: 65
:widths: auto
:align: center

| Column 1     | Column 2      |
| :----------- | :------------ |
| **Text**     |  More Text    |
:::

:::{dropdown} Show code
```:::{table} Table 1: Example table```\
```:width: 65```\
```:widths: auto```\
```:align: center```\
\
```| Column 1     | Column 2      |```\
```| :----------- | :------------ |```\
```| **Text**     |  More Text    |```\
```:::```
:::


### Figures

Here is an example to include a figure, or better two figures, one for dark mode (`:class: only-dark`) and one for light mode (`:class: only-light`).

:::{image} ../../_static/ICON_logo_black.svg
:class: only-light
:height: 100
:width: 200
:::

:::{image} ../../_static/ICON_logo_white.svg
:class: only-dark
:height: 100
:width: 200
:::
*Figure 1: Example figure*

:::{dropdown} Show code
```:::{image} ../../_static/ICON_logo_black.svg```\
```:class: only-light```\
```:height: 100```\
```:width: 200```\
```:::```\
\
```:::{image} ../../_static/ICON_logo_white.svg```\
```:class: only-dark```\
```:height: 100```\
```:width: 200```\
```:::```\
```*Figure 1: Example figure*```
:::

### Namelist References

Since the [Namelist Overview](ref_buildrun_nml) is part of [docs.icon-model.org](https://docs.icon-model.org),
namelist parameters can be referenced from anywhere of [docs.icon-model.org](https://docs.icon-model.org), e.g.

[iforcing](run_nml-iforcing)

:::{dropdown} Show code
```[iforcing](run_nml-iforcing)```
:::

### Admonitions

To highlight a short piece of information, the `admonition-icontheme` can be used

:::{admonition} Some Highlight
:class: admonition-icontheme
Description of the highlight.
:::

:::{dropdown} Show code
```:::{admonition} Some Highlight```

```:class: admonition-icontheme```

```Description of the highlight.```

```:::```
:::

### Highlighting with material icons

In-line highlighting, e.g. for warnings, can be achieved using [material icons](https://fonts.google.com/icons?icon.set=Material+Icons).
See also the [sphinx-design documentation](https://sphinx-design.readthedocs.io/en/latest/badges_buttons.html#material-design-icons) for further information.

{material-regular}`settings;1em;pst-color-secondary` _Some setting_

{material-outlined}`info;2em;pst-color-primary` _Some information_

{material-regular}`warning;2em;pst-color-secondary` _Some warning_

:::{dropdown} Show code
```{material-regular}`settings;1em;pst-color-secondary` _Some settings_```

```{material-outlined}`info;2em;pst-color-primary` _Some information_```

```{material-regular}`warning;2em;pst-color-secondary` _Some warning_```
:::

### Equations

Equations are rendered via JavaScript using [mathjax for sphinx](https://www.sphinx-doc.org/en/master/usage/extensions/math.html#module-sphinx.ext.mathjax).
This allows for the usage of LaTeX syntax to generate equations:

:::{math}
:label: eqtestlabel
\alpha^2 = \frac{\beta}{\Gamma}
:::

Inline math is possible, e.g. {math}`\alpha=42`.

If a label is provided, an equation can be referenced at later point, e.g. equation {eq}`eqtestlabel`.

:::{dropdown} Show code
```:::{math}```\
```:label: eqtestlabel```\
```\alpha^2 = \frac{\beta}{\Gamma}```\
```:::```

```Inline math is possible, e.g. {math}`\alpha=42`.```

```If a label is provided, an equation can be referenced at later point, e.g. equation {eq}`eqtestlabel`.```
:::

### Listing

Apart from the usual hyphon or astrisk syntax in markdown, it is possible to create definition lists which include descriptions of the list items in the following way:

List Item
: Description of the list item.

:::{dropdown} Show code
```List Item```

```: Description of the list item.```
:::

### Intellectual Property Rights

Authors may **optionally** put their name to a section they provided to the documentation. This has no affect on the license. This may serve several purposes:

- Marking the intellectual property rights, thus easing a citation of the respective section
- Provide a contact point for users in case of questions
- A section may not be changed without the consent of the section author

if an author decides to mark the intellectual property right of a section, a macro must be used directly below the section heading for a unified appearance of the web site:

{{author}}A. Name{{endauthor}}

:::{dropdown} Show code
```{{author}}A. Name{{endauthor}}```
:::

Or optionally with institution XYZ:

{{author}}A. Name, XYZ{{endauthor}}

:::{dropdown} Show code
```{{author}}A. Name, XYZ{{endauthor}}```
:::

(ref_docs_guidelines_tor)=
## Terms of Reference

**The task of the docs subsite is to collect:**

- Model documentation (not to replace scientific documentation!)
- Material to run particular configurations
- Introductory material (e.g. Tutorials, Training, ...)

**We do not want:**

- Project websites for particular developments
- Replace or repeat scientific documentation which is published and available.

(ref_docs_guidelines_links)=
## Useful Links

- [**Sphinx**:](https://www.sphinx-doc.org/en/master/) _Read RST or Markdown to generate HTML website_
- [**PyData Sphinx Theme:**](https://pydata-sphinx-theme.readthedocs.io/en/stable/) _A clean, Bootstrap-based Sphinx theme_
- [**MyST:**](https://myst-parser.readthedocs.io/en/) _Markedly Structured Text - Parser_
- [**sphinx{design}:**](https://sphinx-design.readthedocs.io/en/pydata-theme/) _Extension to design screen-reponsive web components_
