/*
A KBase module: samuelbrinkerfunctional_enrichment_update
*/

module samuelbrinkerfunctional_enrichment_update {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_samuelbrinkerfunctional_enrichment_update(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
