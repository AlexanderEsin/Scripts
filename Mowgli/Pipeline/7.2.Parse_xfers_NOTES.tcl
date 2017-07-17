currently in 7.1.Parse_xfers.tcl we are setting a single scenario for output
	however there might be a possibility that two scenarious occur - e.g. both a 3 and a 1.
	need to have some sort of redundant system for associating every event with an individual scenario for downstream processing.