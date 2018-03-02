

Code review
===========


Master and release branches are protected
-----------------------------------------

This means that nobody can push to them directly, all contributions have to be submitted
using merge requests (see below).


How to submit a merge request
-----------------------------

* Before you do anything, create a new branch from the branch that you wish to modify.
* Implement your change and commit to the new branch.
* Push your new branch.
* Go to https://gitlab.com/dirac/dirac and click on the blue button "Create merge request" (top right).
* Describe the change unless it's evident from the commit message.
* Check the box next to "Remove source branch when merge request is accepted".
* Make sure it is directed to the right branch.
* Click green button "Submit merge request".


How to ask for a code review
----------------------------

* Either assign to somebody.
* Or mention somebody via their @username in the discussion. Then they get notified.
* Or let the merge request sit and somebody will pick it up.


How to do code review
---------------------

* Be polite.
* Check that code is submitted to the right branch.
* Check that commit messages are clear.
* Changelog entry needed?
* Documentation needed?
* Tests preserved and new functionality tested?


Further reading
---------------

* https://docs.gitlab.com/ee/gitlab-basics/add-merge-request.html
* https://about.gitlab.com/2017/03/17/demo-mastering-code-review-with-gitlab/
