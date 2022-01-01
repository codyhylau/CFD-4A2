CCN: 5726A
HOW TO RUN EACH OF THE 6 TEST CASES:

1.  Copy the test case geom and flow into files 'geom' and 'flow'. 
    e.g. for test case 0 in the terminal:
    
    cp test0_geom geom
    cp test0_flow flow

    For other test cases, change above command's 0 to whichever test case number from (0-5)
    
    
2.  Edit 'flow' to change CFL and SF as desired. By default they are the best CPU time's CFL and SF. Note for high CFL and/or low SF the solution diverges!


3.  For test cases 4 and 5:
    If geom is copied from the existing test4_geom or test5_geom correctly in step 1, by default the correct flow_guess will be selected in euler.f90. (no need to do anything!)
    Otherwise, check the title in geom is 'supersonic wedge with M1 = 1.6 ' for test case 4 and ' Mach 3 Corner ' for test case 5.
    These are checked in euler.f90, lines 56-68. 'flow_guess_sup is called' for case 4 and 'new_guess' for case 5. The variable 'sup' is set to 1 for case 4.
    
    
4.  The 4 extensions can be enabled/disabled individually in euler.f90 lines 16-19. Setting = 1 to enable, = 0 to disable. By default all 4 are enabled.


5.  To compile, just run it normally. i.e. in the terminal:

    make clean
    make
    ./Euler
    
    