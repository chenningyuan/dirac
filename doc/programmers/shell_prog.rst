Shell programming
=================

It is hard to find good references on true bourne shell programming. An
extensive list of the differences between bash and sh can be found
`here <http://www.faqs.org/faqs/unix-faq/shell/bash>`_.

Nevertheless, we are working with python language <www.python.org>, and all new scripts
are supposed to be written in this language.

When writing python scripts, please bear in mind that these have to be platform general,
i.e. working properly under both Windows and Unix systems. Please avoid using specific Unix commands
in subprocesses, as these shall be usable only under Unix family of operating systems.
This wish holds specially for *setup*, *runtest* and *pam* scripts.


