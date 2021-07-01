Review Checklist:

- Code review is for code improvement, not a merge/not merge decision
- Be thankful for code contributions
- Make sure to maintain code quality, merge is a privilege, not an automatism

Start with the large scale picture:
===================================

1. Is the change of possible interest for other users?

This should by default be answered with 'yes'. If one person needs it it is likely useful for others as well. However, there are changes that are too specific, or too disruptive to be considered for integration (say changing ASPECT into a finite-difference code not based on deal.II, or replacing deal.II with a different FE library).

2. Does the current structure fit into ASPECT's design?

Is the plugin in the right place? Does it belong into this plugin system or directory? Keep in mind that ASPECT's structure can be adaped if necessary. Large-scale reorganization is painful, but preferrable to integrating hacks that deteriorate the quality of the code base.

3. Is backward compatibility preserved?

Breaking backward compatibility is acceptable if absolutely necessary and for a good and sufficiently important reason, but it should be avoided otherwise.

4. Are examples in the right folder?

We separate examples into benchmarks (accuracy/convergence checks), cookbooks (instructional models), and tests (regression tests with a known 'good' solution).


Then consider small-scale questions:
====================================

- Does the code look like ASPECT code? Does it follow coding and naming conventions?
- Are indentation, comments, documentation structured correctly?
- Is encapsulation preserved?
- Is the purpose clear to users?
- Is it possible to significantly improve the code with reasonable effort?
- Does the PR have at least a test?
- Does the PR have a changelog entry (for changes that are significant for users)?
