:orphan:
 

How to add/move/remove sources
------------------------------

The ``pamadm`` script automatically
assembles the list of sources.
Since ``pamadm`` is run every time you hit
"make", you normally don't need to do anything
special after adding or moving or removing
a file - it should be automatically detected
and "just work" (TM).


How to avoid that sources get added
-----------------------------------

Add the file(s) to ``files_never_compiled`` in ``pamadm``.


How to avoid that an entire directory gets added
------------------------------------------------

Add the directory to ``directories_never_compiled`` in ``pamadm``.
Alternatively create a ``pamadm`` module. With this you
can toggle whether a directory gets built or not.
