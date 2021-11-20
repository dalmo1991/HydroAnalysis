.. _contribute:

Software organization and contribution
======================================

The HydroAnalysis framework comprises the following components:

- **Source code**: Latest version of all the code necessary to use the
  framework. The source code would normally be accessed only by advanced
  users, e.g. to understand the internal organization of the framework, to install
  manually the latest version, to extend the framework with new
  functionality, etc.
- **Packaged release**: Latest stable version of the framework available for
  users.
- **Documentation**: Detailed explanation of the framework.

The source code, documentation, and examples are part of the official repository
of HydroAnalysis hosted on `GitHub <https://github.com/dalmo1991/HydroAnalysis>`_.
A user who wishes to read the source code and/or modify any aspect of
HydroAnalysis (source code, documentation, and examples) can do it using GitHub.

New releases of the software are available from the official Python Package
Index (PyPI), where HydroAnalysis has a
`dedicated page <https://pypi.org/project/HydroAnalysis/>`_.

The documentation builds automatically from the
`source folder <https://github.com/dalmo1991/HydroAnalysis/tree/master/doc>`_ on
GitHub and is published online in
`Read the Docs <https://HydroAnalysis.readthedocs.io/>`_.

Contributions
-------------

Contributions to the framework can be made in the following ways:

- Submit issues on bugs, desired features, etc;
- Solve open issues;
- Extend the documentation with new demos and examples;
- Extend and/or modify the framework;
- Use and cite the framework in your publications.

Code contribution by external users will be mainly additive (i.e., adding new
components, as illustrated in :ref:`build_element` and :ref:`customize_components`)
and should include also appropriate testing (:ref:`tests`).

Contributors will maintain authorship of the contributed code and are invited
to include, in all files, their contact information to facilitate future
collaboration. The authors and maintainers of HydroAnalysis will undertake a basic
inspection of the contributed code to identify any quality issues.

The typical workflow that should be followed when contributing to a GitHub
project is described
`here <https://www.dataschool.io/how-to-contribute-on-github/>`_.

In summary,
the following steps should be followed:

1. Fork the HydroAnalysis repository to the user GitHub account;
2. Clone the fork on the user computer;
3. Modify the code, commit the changes, and push them to the GitHub fork of
   HydroAnalysis;
4. Make a pull request on GitHub to the HydroAnalysis repository.

Branching scheme of the GitHub repository
.........................................

Updates to HydroAnalysis are made directly in the branch :code:`master`, which
is the most up-to-date branch. The branch :code:`release` is used only
for the staging of new software releases and, therefore, code should not be
pushed directly to it.

When a code update is merged from :code:`master` to :code:`release`, a
new version of the package is automatically released on PyPI. Remember to update
the version number in the :code:`setup.py` file to avoid conflicts.

Developers are free to create new branches, but pull requests must be directed to
:code:`master` and not to :code:`release`.

Documentation and examples are generated from the :code:`master`
branch.
