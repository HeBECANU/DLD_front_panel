This is a upgraded version of the Classic dld_front_panel.
-suports reading in TXY files
	-checks that they were modified after the TDC data file. IF not reconvert
-faster conversion
	-thanks to not using shit code (we dont have plate pulses anyway)
-rewritten monitor code
 	-seemles free run/ normal operation
	- ability to wait for TXY files (for use with auto_convert on TDC computer)
	-more text output for monitoring
-pressign enter in the TOF min/max bin boxes will refresh the TOF window
-enter in filepath and file number will read data
-reading data will auto update TOF and 2d
-rewriten 2d code
    -2d plot can show in XT and YT mode
    -all 3 combined plot
    -no aliasing of bins
    -way faster using internal hist3
-rewritten TOF code
    -faster using internal hist
    -TOF reads in kHz
-3d plot
    -allows rotation


