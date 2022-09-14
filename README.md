# RAREdar
Cleaned up version of the newRAREdar

# TODOs

- update this README with brief instructions about how to run. Rename the screenshots and include them in the README.
- Create directories to organize your files: output/ input/ src/
- Merge this branch into `main` when you're done. [Done]

Comments/questions for `RAREdar.py`:
- Is `ReverseList` the reverse or the reverse complement of `MotifList`?  Add comments that carefully explains this. [Done]
- Consider renaming `Gene_Sequence_list.fa` with the species name, in case we ever want to run it on another file. [Done]
- Speedup suggestion: in your FOR loop, first check if `window==repeat` and THEN check the DRList/RDRList, etc. If the window `!=` repeat then you can skip that whole FOR loop.
- Your `if/elif` structure is good if you want to count the number of CODING RARE, COMPLEMENTARY RARAE, REVERSED RARE, or REVERSED COMPLEMENTARY RARE hits.  If you don't care about this, though, that whole structure could be turned into a single `if` with the logic connected by `or`s.  
