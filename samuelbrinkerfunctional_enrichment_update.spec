/*
A KBase module: samuelbrinkerfunctional_enrichment_update
*/

module samuelbrinkerfunctional_enrichment_update {
    /* A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int boolean;

    /* An X/Y/Z style reference
    */
    typedef string obj_ref;

    /*
      required params:
      feature_set_ref: FeatureSet object reference
      workspace_name: the name of the workspace it gets saved to

      optional params:
      propagation: includes is_a relationship to all go terms (default is 1)
      filter_ref_features: filter reference genome features with no go terms (default is 0)
      statistical_significance: parameter for statistical significance. Select one from left_tailed, right_tailed or two_tailed (default is left_tailed)
      ignore_go_term_not_in_feature_set: ignore Go term analysis if term is not associated with FeatureSet (default is 1)
    */
    typedef structure {
        obj_ref feature_set_ref;
        string workspace_name;
        boolean propagation;
        boolean filter_ref_features;
        string  statistical_significance;
        boolean ignore_go_term_not_in_feature_set;
        string orthology_type;

    } AppInput;
    typedef structure{
        string result_directory;
        string report_name;
        string report_ref;  }ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_samuelbrinkerfunctional_enrichment_update(AppInput params) returns (ReportResults output) authentication required;

};
