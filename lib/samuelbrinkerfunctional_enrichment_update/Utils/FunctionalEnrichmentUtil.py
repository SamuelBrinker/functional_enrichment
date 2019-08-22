import csv
import errno
import json
import os
import re
import time
import uuid
import zipfile
import fisher

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

#from installed_clients.DataFileUtilClient import DataFileUtil
#from installed_clients.GenomeSearchUtilClient import GenomeSearchUtil
#from installed_clients.KBaseReportClient import KBaseReport
#from GenomeSearchUtil.GenomeSearchUtilImpl import GenomeSearchUtil
#from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.GenomeSearchUtilClient import GenomeSearchUtil
#from DataFileUtil.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
#from Workspace.WorkspaceClient import Workspace as Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class FunctionalEnrichmentUtil:

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _validate_run_fe1_params(self, params):
        """
        _validate_run_fe1_params:
                validates params passed to run_fe1 method
        """

        log('start validating run_fe1 params')

        # check for required parameters
        for p in ['feature_set_ref', 'workspace_name']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _generate_report(self, enrichment_map, result_directory, workspace_name,
                         feature_id_go_id_list_map, feature_set_ids, genome_ref,
                         go_id_parent_ids_map, feature_ids):
        """
        _generate_report: generate summary report
        """

        log('start creating report')

        output_files = self._generate_output_file_list(result_directory,
                                                       enrichment_map,
                                                       feature_id_go_id_list_map,
                                                       feature_set_ids,
                                                       genome_ref,
                                                       go_id_parent_ids_map,
                                                       feature_ids)

        output_html_files = self._generate_html_report(result_directory,
                                                       enrichment_map)

        report_object_name = 'functional_enrichment_update_report_' + str(uuid.uuid4())
        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 333,
                         'report_object_name': report_object_name}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_supporting_files(self, result_directory, enrichment_map,
                                   feature_id_go_id_list_map, feature_set_ids, genome_ref,
                                   go_id_parent_ids_map, feature_ids):
        """
        _generate_supporting_files: generate varies debug files
        """
        supporting_files = list()

        feature_id_go_ids_map_file = os.path.join(result_directory, 'feature_id_go_ids_map.txt')
        go_id_genome_feature_ids_map_file = os.path.join(result_directory,
                                                         'go_id_genome_feature_ids_map.txt')
        go_id_set_feature_ids_map_file = os.path.join(result_directory,
                                                      'go_id_feature_set_feature_ids_map.txt')
        feature_ids_file = os.path.join(result_directory, 'feature_ids.txt')
        feature_set_ids_file = os.path.join(result_directory, 'feature_set_ids.txt')
        fisher_variables_file = os.path.join(result_directory, 'fisher_variables.txt')
        genome_info_file = os.path.join(result_directory, 'genome_info.txt')
        go_id_parent_ids_map_file = os.path.join(result_directory, 'go_id_parent_ids_map.txt')

        supporting_files.append(feature_id_go_ids_map_file)
        supporting_files.append(go_id_genome_feature_ids_map_file)
        supporting_files.append(feature_ids_file)
        supporting_files.append(feature_set_ids_file)
        supporting_files.append(fisher_variables_file)
        supporting_files.append(genome_info_file)
        supporting_files.append(go_id_parent_ids_map_file)
        supporting_files.append(go_id_set_feature_ids_map_file)

        total_feature_ids = list(feature_id_go_id_list_map.keys())
        feature_ids_with_feature = []
        for feature_id, go_ids in feature_id_go_id_list_map.items():
            if isinstance(go_ids, list):
                feature_ids_with_feature.append(feature_id)
        genome_name = self.ws.get_object_info3({'objects':
                                                [{'ref': genome_ref}]})['infos'][0][1]

        with open(go_id_parent_ids_map_file, 'w') as go_id_parent_ids_map_file:
            for go_id, parent_ids in go_id_parent_ids_map.items():
                go_id_parent_ids_map_file.write(f'{go_id}: {", ".join(parent_ids)}\n')

        with open(genome_info_file, 'w') as genome_info_file:
            genome_info_file.write(f'genome_name: {genome_name}\n')
            genome_info_file.write(f'features: {len(total_feature_ids)}\n')
            genome_info_file.write(f'features with term: {len(feature_ids_with_feature)}')

        with open(feature_set_ids_file, 'w') as feature_set_ids_file:
            feature_set_ids_file.write('\n'.join(feature_set_ids))

        with open(feature_id_go_ids_map_file, 'w') as feature_id_go_ids_map_file:
            with open(feature_ids_file, 'w') as feature_ids_file:
                for feature_id, go_ids in feature_id_go_id_list_map.items():
                    feature_ids_file.write(f'{feature_id} {feature_id in feature_set_ids}\n')
                    if isinstance(go_ids, str):
                        feature_id_go_ids_map_file.write(f'{feature_id} {go_ids}\n')
                    else:
                        feature_id_go_ids_map_file.write(f'{feature_id} {", ".join(go_ids)}\n')

        with open(go_id_genome_feature_ids_map_file, 'w') as go_id_genome_feature_ids_map_file:
            with open(go_id_set_feature_ids_map_file, 'w') as go_id_set_feature_ids_map_file:
                with open(fisher_variables_file, 'w') as fisher_variables_file:
                    for go_id, go_info in enrichment_map.items():
                        mapped_features = go_info.get('mapped_features')
                        fs_mapped_features = list(set(mapped_features).intersection(
                          feature_set_ids))
                        mapped_features_line = f'{go_id}: {", ".join(mapped_features)}\n'
                        go_id_genome_feature_ids_map_file.write(mapped_features_line)

                        set_mapped_features_line = f'{go_id}: {", ".join(fs_mapped_features)}\n'
                        go_id_set_feature_ids_map_file.write(set_mapped_features_line)
                        a_value = go_info.get('num_in_subset_feature_set')
                        b_value = len(feature_set_ids) - a_value
                        c_value = len(mapped_features) - a_value
                        d_value = len(feature_ids) - len(feature_set_ids) - c_value
                        p_value = go_info.get('raw_p_value')
                        fisher_variables_file.write(
                            f'{go_id} a:{a_value} b:{b_value} c:{c_value} d:{d_value} ')
                        fisher_variables_file.write(f'p_value:{p_value}\n')

        result_file = os.path.join(result_directory, 'supporting_files.zip')
        with zipfile.ZipFile(result_file, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as zip_file:
            for supporting_file in supporting_files:
                zip_file.write(supporting_file,
                               os.path.basename(supporting_file))

        return [{'path': result_file,
                 'name': os.path.basename(result_file),
                 'label': os.path.basename(result_file),
                 'description': 'GO term functional enrichment supporting files'}]

    def _generate_output_file_list(self, result_directory, enrichment_map,
                                   feature_id_go_id_list_map, feature_set_ids, genome_ref,
                                   go_id_parent_ids_map, feature_ids):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """

        log('start packing result files')
        output_files = list()

        result_file = os.path.join(result_directory, 'functional_enrichment.csv')
        with open(result_file, 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['term_id', 'term', 'ontology', 'num_in_feature_set',
                             'num_in_ref_genome', 'raw_p_value', 'adjusted_p_value'])
            for key, value in enrichment_map.items():
                writer.writerow([key, value['go_term'], value['namespace'],
                                 value['num_in_subset_feature_set'],
                                 value['num_in_ref_genome'], value['raw_p_value'],
                                 value['adjusted_p_value']])

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'GO term functional enrichment'})

        supporting_files = self._generate_supporting_files(result_directory,
                                                           enrichment_map,
                                                           feature_id_go_id_list_map,
                                                           feature_set_ids,
                                                           genome_ref,
                                                           go_id_parent_ids_map,
                                                           feature_ids)
        output_files += supporting_files

        return output_files

    def _generate_html_report(self, result_directory, enrichment_map):
        """
        _generate_html_report: generate html summary report
        """

        log('start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        enrichment_table = ''
        data = csv.DictReader(open(os.path.join(result_directory, 'functional_enrichment.csv')),
                              delimiter=',')
        sortedlist = sorted(data, key=lambda row: (float(row['adjusted_p_value']),
                                                   float(row['raw_p_value']),
                                                   float(row['num_in_ref_genome'])),
                            reverse=False)

        for row in sortedlist:
            # if row['num_in_feature_set'] != '0':
            enrichment_table += f'<tr><td>{row["term_id"]}</td>'
            enrichment_table += f'<td>{row["term"]}</td>'
            enrichment_table += f'<td>{row["ontology"]}</td>'
            enrichment_table += f'<td>{row["num_in_feature_set"]}</td>'
            enrichment_table += f'<td>{row["num_in_ref_genome"]}</td>'
            enrichment_table += f'<td>{float(row["raw_p_value"]):.3g}</td>'
            enrichment_table += f'<td>{float(row["adjusted_p_value"]):.3g}</td></tr>'

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<tr>Enrichment_Table</tr>',
                                                          enrichment_table)
                result_file.write(report_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for Functional Enrichment App'})
        return html_report

    def _get_go_maps_from_genome(self, genome_ref,orthology_type):
        """
        _search_genome: search genome data
        """

        log('start parsing GO terms from genome')

        feature_num = self.gsu.search({'ref': genome_ref})['num_found']

        genome_features = self.gsu.search({'ref': genome_ref,
                                           'limit': feature_num,
                                           'sort_by': [['feature_id', True]]})['features']

        feature_id_go_id_list_map = {}
        go_id_feature_id_list_map = {}
        go_id_go_term_map = {}
        feature_id_feature_info_map = {}

        for genome_feature in genome_features:
            feature_id = genome_feature.get('feature_id')
            feature_func = genome_feature.get('function')
            feature_type = genome_feature.get('feature_type')
            ontology_terms = genome_feature.get('ontology_terms')
            feature_id_feature_info_map.update({feature_id: {'function': feature_func,
                                                             'feature_type': feature_type}})

            go_id_list = []
            if ontology_terms:
                for ontology_id, ontology_term in ontology_terms.items():
                    if orthology_type=='GO' and re.match('[gG][oO]\:.*', ontology_id):
                        go_id_go_term_map.update({ontology_id: ontology_term})
                        go_id_list.append(ontology_id)
                    elif orthology_type=='MetaCyc' and  re.match('RXN.*', ontology_id):
                        go_id_go_term_map.update({ontology_id: ontology_term})
                        go_id_list.append(ontology_id)
                    elif orthology_type=='Kegg' and  re.match('R\d.*', ontology_id):
                        go_id_go_term_map.update({ontology_id: ontology_term})
                        go_id_list.append(ontology_id)
                    elif orthology_type=='EC' and re.match('\d\..*', ontology_id) and not re.match('RXN.*', ontology_id):
                        go_id_go_term_map.update({ontology_id: ontology_term})
                        go_id_list.append(ontology_id)
                    elif orthology_type=='all_terms' and (re.match('[gG][oO]\:.*', ontology_id) or
                                    re.match('RXN.*', ontology_id) or re.match('R\d.*', ontology_id)
                                    or re.match('\d\..*', ontology_id)):
                        go_id_go_term_map.update({ontology_id: ontology_term})
                        go_id_list.append(ontology_id)

            if go_id_list:
                feature_id_go_id_list_map.update({feature_id: go_id_list})

                for go_id in go_id_list:
                    if go_id in go_id_feature_id_list_map:
                        feature_ids = go_id_feature_id_list_map.get(go_id)
                        feature_ids.append(feature_id)
                        go_id_feature_id_list_map.update({go_id: feature_ids})
                    else:
                        go_id_feature_id_list_map.update({go_id: [feature_id]})
            else:
                feature_id_go_id_list_map.update({feature_id: 'Unlabeled'})

        return (feature_id_go_id_list_map, go_id_feature_id_list_map,
                go_id_go_term_map, feature_id_feature_info_map)

    def _process_feature_set(self, feature_set_ref):
        """
        _process_feature_set: process FeatureSet object

        return:
        genome_ref: reference Genome object ref
        feature_set_ids: FeatureSet feature ids
        """

        log('start processing FeatureSet object')

        feature_set_data = self.ws.get_objects2({'objects':
                                                 [{'ref': feature_set_ref}]})['data'][0]['data']
        feature_elements = feature_set_data['elements']
        feature_set_ids = []
        genome_ref_array = []
        for feature_id, genome_refs in feature_elements.items():
            feature_set_ids.append(feature_id)
            genome_ref_array += genome_refs

        if len(set(genome_ref_array)) > 1:
            error_msg = 'FeatureSet has multiple reference Genomes: {}'.format(genome_ref_array)
            raise ValueError(error_msg)

        return feature_set_ids, genome_ref_array[0]


        #####################################################
    def _get_immediate_parents(self, ontology_hash, go_id, is_a_relationship=True,
                               regulates_relationship=True, part_of_relationship=False):
        """
        _get_immediate_parents: get immediate parents go_ids for a given go_id
        """
        parent_ids = []
        antology_info = ontology_hash.get(go_id, {})

        if is_a_relationship:
            is_a_parents = antology_info.get('is_a')
            if is_a_parents:
                for parent_string in is_a_parents:
                    is_a_parent_id = parent_string.split('!')[0][:-1]
                    parent_ids.append(is_a_parent_id)

        if regulates_relationship:
            relationship = antology_info.get('relationship')
            if relationship:
                for relationship_string in relationship:
                    if relationship_string.split(' ')[0] == 'regulates':
                        parent_ids.append(relationship_string.split(' ')[1])

        if part_of_relationship:
            relationship = antology_info.get('relationship')
            if relationship:
                for relationship_string in relationship:
                    if relationship_string.split(' ')[0] == 'part_of':
                        parent_ids.append(relationship_string.split(' ')[1])

        return parent_ids

    def _get_kegg_parents(self, kegg):
        ### This file is a hierarchical map of kegg IDs
        with open('/kb/module/data/br08201.json', 'r') as f:
            kegg_reactions = json.load(f)
            f.close()

        #kegg=list(ontology_hash.keys())
        print("Getting Kegg parents")
        kegg_id_reaction_ids={}
         ####change
        kegg_id_reaction_ids={}
         ####change

        for id in kegg:
            true_reactions=[]
            parrent_ids=[]
            base_level=''
            early_break=False
            while early_break==False:
                for map in kegg_reactions['children']:
                    for children_1 in map['children']:
                        for children_2 in children_1['children']:
                            for children_3 in children_2['children']:
                                if 'children' in children_3:
                                    for children_4 in children_3['children']:
                                        if id in children_4['name']:
                                            base_level=children_3['name']
                                            early_break=True
                early_break=True
            if base_level!='':
                split_base=base_level.split('.')
                while len(split_base)<4:
                    split_base.append('-')
                i=len(split_base)-1
                while i >=0:
                    if '-' not in split_base[i]:
                        split_base[i]='-'
                        possible_parent='.'.join(split_base)
                        for map in kegg_reactions['children']:
                            for children_1 in map['children']:
                                for children_2 in children_1['children']:
                                    for children_3 in children_2['children']:
                                        if 'children' in children_3 and possible_parent in children_3['name']:
                                            for children_4 in children_3['children']:
                                                parrent_ids.append(children_4['name'])
                    i-=1
                for other_reaction in kegg:
                    for parent in parrent_ids:
                        if other_reaction in parent and other_reaction not in true_reactions:
                            true_reactions.append(other_reaction)
                kegg_id_reaction_ids[id]=true_reactions

        with open('/kb/module/data/br08901.json', 'r') as f:
            kegg_reactions = json.load(f)
            f.close()
        for id in kegg:
            true_reactions=[]
            parrent_ids=[]
            base_level=''
            early_break=False
            for cat in  kegg_reactions['children']:
                for cat2 in cat['children']:
                    for cat3 in cat2['children']:
                        if id.replace('R','') in str(cat3['name']):
                            for cat_check in cat2['children']:
                                for all_id in kegg:
                                    if all_id.replace('R','') in str(cat_check['name']) and all_id !=id:
                                        if id not in kegg_id_reaction_ids.keys():
                                            kegg_id_reaction_ids[id]=[all_id]
                                        elif all_id not in kegg_id_reaction_ids[id]:
                                            kegg_id_reaction_ids[id]=kegg_id_reaction_ids[id] +[all_id]

        return kegg_id_reaction_ids

    def _get_ec_parents(self, ecs):

        #ecs=list(ontology_hash.keys())
        ec_id_parent_id={}
        print("Getting EC Parents")
        for ec in ecs:
            split_ec=ec.split('.')
            while len(split_ec)<4:
                split_ec.append('-')
            i=len(split_ec)-1
            parents=[]
            while i >=0:
                if '-' not in split_ec[i]:
                    split_ec[i]='-'
                    possible_parent='.'.join(split_ec)
                    if possible_parent in ecs and possible_parent not in parents: ####### Change if part of enrichment
                        parents.append(possible_parent)
                i-=1
            ec_id_parent_id[ec]=parents

        return ec_id_parent_id


    def _get_metacyc_reactions(self, metacyc):

        def pathways_to_reactions(self, pathway):
            pathway_reacionts=[]
            for chunk in metacyc_raw.split('</Reaction>'):
                if pathway in chunk:
                    pathway_reacionts.append(chunk.split("<Reaction ID='")[1].split("' orgid='")[0])
            return pathway_reacionts


        with open('/kb/module/data/metacyc_reactions.xml','r') as file:
            metacyc_raw=file.read()
            file.close()

        #metacyc=list(ontology_hash.keys())
        metacyc_id_pathway_id={}
        print("Getting MetaCyc parents")
        for meta in metacyc:
            metacyc_id_pathway_id[meta]=[]
            pathway_reacionts=[]
            pathways=[]
            #Cycle through each reaction until metacyc ID is found
            for chunk in metacyc_raw.split('</Reaction>'):
                if "<Reaction ID='"+meta+"'" in chunk:
                    #Extract names of pathways it is present in
                    if "<in-pathway>" in chunk:
                        split_pathways=chunk.split('<in-pathway>')[1].split("</in-pathway>")[0].split('\n')
                        for line in split_pathways:
                            if 'Pathway' in line:
                                pathways.append(line.split("<Pathway resource='getxml?")[1].split("' orgid")[0])
            for path in pathways:
                pathway_reacionts=self.pathways_to_reactions(path)
                for reaction in pathway_reacionts:
                    if reaction in metacyc and reaction !=meta:
                        metacyc_id_pathway_id[meta].append(reaction)

        return metacyc_id_pathway_id


    def _parents_from_all_terms(self, ontology_hash, is_a_relationship=True, regulates_relationship=True, part_of_relationship=False):

        all_ids=list(ontology_hash.keys())
        go_ids=[]
        kegg_ids=[]
        metacyc_ids=[]
        ec_ids=[]

        all_parents=[]
        all_reactions=[]
        go_reactions=[]
        kegg_reactions=[]
        ec_reactions=[]
        metacyc_reactions=[]

        go_parents={}
        kegg_parents={}
        ec_parents={}
        metacyc_parents={}

        ids_with_parents={}
        print("Getting All Parents")
        for id in all_ids:
            if 'GO:' in id.upper():
                go_ids.append(id)
            elif 'RXN' in id.upper():
                metacyc_ids.append(id)
            elif re.match('R\d.*', id):
                kegg_ids.append(id)
            elif re.match('\d\..*', id):
                ec_ids.append(id)

        if go_ids!=[]:
            go_parents={}
            print("Getting GO Parents")
            for go_id in go_ids:
                go_parents[go_id]=self._get_go_parents(ontology_hash, go_id,
                                                         is_a_relationship, regulates_relationship,
                                                         part_of_relationship)
            with open('/kb/module/data/GO.ModelSEED.json', 'r') as f:
                table = json.load(f)
                f.close()
            print("Translating GO terms")
            go_reactions= self._translate_terms( list(go_parents.keys()) , table )

        if kegg_ids!=[]:
            kegg_parents= self._get_kegg_parents(kegg_ids)
            with open('/kb/module/data/KEGG_RXN.ModelSEED.json', 'r') as f:
                table = json.load(f)
                f.close()
            print("Translating Kegg terms")
            kegg_reactions= self._translate_terms( list(kegg_parents.keys()) , table )

        if ec_ids!=[]:
            ec_parents= self._get_ec_parents(ec_ids)
            with open('/kb/module/data/EBI_EC.ModelSEED.json', 'r') as f:
                table = json.load(f)
                f.close()
            print("Translating EC terms")
            ec_reactions= self._translate_terms( list(ec_parents.keys()) , table )

        if metacyc_ids!=[]:
            metacyc_parents= self._get_metacyc_reactions(metacyc_ids)
            with open('/kb/module/data/Metacyc_RXN.ModelSEED.json', 'r') as f:
                table = json.load(f)
                f.close()
            print("Translating MetaCyc terms")
            metacyc_reactions= self._translate_terms( list(metacyc_parents.keys()) , table )

        all_parents.append(go_parents)
        all_parents.append(kegg_parents)
        all_parents.append(ec_parents)
        all_parents.append(metacyc_parents)

        all_reactions.append(go_reactions)
        all_reactions.append(kegg_reactions)
        all_reactions.append(ec_reactions)
        all_reactions.append(metacyc_reactions)

        print("Comparing reactions")
        for bunch in all_reactions:
            bunch_term= all_reactions.index(bunch)
            for comparison in all_reactions:
                comparison_term= all_reactions.index(comparison)
                if bunch!=comparison:
                    for bunch_reaction in bunch:
                        for comparison_reaction in comparison:
                            if bunch_reaction==comparison_reaction and bunch_reaction!='None' and comparison_reaction!='':

                                bunch_term_position=bunch.index(bunch_reaction)
                                comparison_term_position=comparison.index(comparison_reaction)

                                bunch_id=all_parents[bunch_term].keys()[bunch_term_position]
                                comparison_id=all_parents[comparison_term].keys()[comparison_term_position]

                                all_parents[bunch_term][bunch_id]=all_parents[bunch_term][bunch_id]+all_parents[comparison_term][comparison_id]
                                all_parents[bunch_term][bunch_id].append(comparison_id)

        for bunch in all_parents:
            ids_with_parents.update(bunch)


        return ids_with_parents


    def _translate_terms(self, ids, table):
        translated_ids=[]
        for id in ids:
            full_reaction=[]
            if id.replace('GO:','').replace('META:','') in list(table['translation'].keys()):
                for reaction in table['translation'][id.replace('GO:','').replace('META:','')]['equiv_terms']:
                    if str(reaction['equiv_term']) not in full_reaction:
                        full_reaction.append(str(reaction['equiv_term']))
            translated_ids.append(' '.join(full_reaction))
        return translated_ids




    def _get_go_parents(self, ontology_hash, go_id, is_a_relationship=True,
                            regulates_relationship=True, part_of_relationship=False):

        grand_parent_ids=[]
        parent_ids = self._get_immediate_parents(ontology_hash, go_id,
                                                 is_a_relationship, regulates_relationship,
                                                 part_of_relationship)
        if parent_ids:
            grand_parent_ids = parent_ids
            for parent_id in parent_ids:
                grand_parent_ids += self._get_go_parents(ontology_hash, parent_id,
                                                         is_a_relationship, regulates_relationship,
                                                         part_of_relationship)

        #fetch_result[go_id]=list(set(grand_parent_ids))
        return list(set(grand_parent_ids))


    def _fetch_all_parents_go_ids(self, ontology_hash, go_ids,orthology_type, is_a_relationship=True,
                                  regulates_relationship=True, part_of_relationship=False):
        '''
        _fetch_all_parents_go_ids: recusively fetch all parent go_ids
        '''
        fetch_result={}
        if orthology_type=='GO':
            print("Getting GO Parents")
            for go_id in go_ids:
                fetch_result[go_id]=self._get_go_parents(ontology_hash, go_id, is_a_relationship, regulates_relationship, part_of_relationship)

            return fetch_result

        elif orthology_type=='Kegg':
            kegg=list(ontology_hash.keys())
            return self._get_kegg_parents(go_ids)
        elif orthology_type=='EC':
            #ecs=list(ontology_hash.keys())
            return self._get_ec_parents(go_ids)
        elif orthology_type=='MetaCyc':
            #metacyc=list(ontology_hash.keys())
            return self._get_metacyc_reactions(go_ids)
        elif orthology_type=='all_terms':
            return self._parents_from_all_terms(ontology_hash,  is_a_relationship,
                                          regulates_relationship, part_of_relationship)
        else:
            return {go_id: []}

    def _generate_parent_child_map(self, ontology_hash, go_ids,orthology_type,
                                   is_a_relationship=True,
                                   regulates_relationship=True,
                                   part_of_relationship=False):
        """
        _generate_parent_child_map: fetch parent go_ids for given go_id
        """

        log('start fetching parent go_ids')
        start = time.time()

        go_id_parent_ids_map = {}

        fetch_result = self._fetch_all_parents_go_ids(ontology_hash, go_ids,orthology_type, is_a_relationship, regulates_relationship, part_of_relationship)

        go_id_parent_ids_map.update(fetch_result)

        end = time.time()
        print(f'used {end - start:.2f} s')

        return go_id_parent_ids_map

    def _round(self, number, digits=3):
        """
        round number to given digits
        """

        round_number = format(number, f'.{digits}g')

        return round_number

    def __init__(self, config):
        self.ws_url = config['workspace-url']
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.scratch = config['scratch']
        self.dfu = DataFileUtil(self.callback_url)
        self.gsu = GenomeSearchUtil(self.callback_url)
        self.ws = Workspace(self.ws_url, token=self.token)

    def run_fe1(self, params):
        """
        run_fe1: Functional Enrichment One

        required params:
        feature_set_ref: FeatureSet object reference
        workspace_name: the name of the workspace it gets saved to

        optional params:
        propagation: includes is_a relationship to all go terms (default is 1)
        filter_ref_features: filter reference genome features with no go terms (default is 0)
        statistical_significance: parameter for statistical significance.
                                  Select one from left_tailed, right_tailed or two_tailed
                                  (default is left_tailed)
        ignore_go_term_not_in_feature_set: ignore Go term analysis if term is not associated with
                                           FeatureSet (default is 1)

        return:
        result_directory: folder path that holds all files generated by run_deseq2_app
        report_name: report name generated by KBaseReport
        report_ref: report reference generated by KBaseReport
        """
        log('--->\nrunning FunctionalEnrichmentUtil.run_fe1\n' +
            f'params:\n{json.dumps(params, indent=1)}')

        self._validate_run_fe1_params(params)
        print(params)

        propagation = params.get('propagation', True)
        orthology_type = params.get('orthology_type', 'GO') ######################################################## Needs to be added
        filter_ref_features = params.get('filter_ref_features', False)
        statistical_significance = params.get('statistical_significance', 'left_tailed')
        ignore_go_term_not_in_feature_set = params.get('ignore_go_term_not_in_feature_set', True)

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)


        print(orthology_type)

        feature_set_ids, genome_ref = self._process_feature_set(params.get('feature_set_ref'))

        ############################################
        # Assuming Kegg / EC data would be stored the
        # same way as GO data for below
        #############################################


        (feature_id_go_id_list_map,
         go_id_feature_id_list_map,
         go_id_go_term_map,
         feature_id_feature_info_map) = self._get_go_maps_from_genome(genome_ref,orthology_type)

        if not len(feature_id_go_id_list_map):
            raise ValueError("No features in the referenced genome ({}) contain ontology mappings"
                             .format(genome_ref))

        unknown_feature_ids = set(feature_set_ids) - set(feature_id_feature_info_map.keys())
        if unknown_feature_ids:
            raise ValueError("The specified feature set contains {} feature ids which are not "
                             "present referenced genome".format(genome_ref))

        if filter_ref_features:
            log('start filtering features with no term')
            feature_ids = []
            for feature_id, go_ids in feature_id_go_id_list_map.items():
                if isinstance(go_ids, list):
                    feature_ids.append(feature_id)
        else:
            feature_ids = list(feature_id_go_id_list_map.keys())

        ontology_hash = dict()


        ############################################################
        #As far as I know, KBaseOntology's gene_ontology is not working
        ############################################################


        ontologies = self.ws.get_objects([{'workspace': 'KBaseOntology',
                                           'name': 'gene_ontology'},
                                          {'workspace': 'KBaseOntology',
                                           'name': 'plant_ontology'}])
        ontology_hash.update(ontologies[0]['data']['term_hash'])
        ontology_hash.update(ontologies[1]['data']['term_hash'])

        if propagation:
            go_id_parent_ids_map = self._generate_parent_child_map(ontology_hash,
                                                                   list(go_id_go_term_map.keys()),orthology_type,
                                                                   regulates_relationship=False)
        else:
            go_id_parent_ids_map = {}
            for go_id in go_id_go_term_map.keys():
                go_id_parent_ids_map.update({go_id: []})

        log('including parents to feature id map')
        for go_id, parent_ids in go_id_parent_ids_map.items():
            mapped_features = go_id_feature_id_list_map.get(go_id)

            for parent_id in parent_ids:
                parent_mapped_features = go_id_feature_id_list_map.get(parent_id)

                if not parent_mapped_features:
                    parent_mapped_features = []

                if mapped_features:
                    parent_mapped_features += mapped_features

                go_id_feature_id_list_map.update({parent_id: list(set(parent_mapped_features))})

        log('start calculating p-values')
        enrichment_map = {}
        go_info_map = {}
        all_raw_p_value = []
        pos = 0
        for go_id, go_term in go_id_go_term_map.items():
            mapped_features = go_id_feature_id_list_map.get(go_id)
            # in feature_set matches go_id
            a = len(set(mapped_features).intersection(feature_set_ids))
            # ignore go term analysis if not associated with FeatureSet
            if ignore_go_term_not_in_feature_set and a == 0:
                continue
            # in feature_set doesn't match go_id
            b = len(feature_set_ids) - a
            # not in feature_set matches go_id
            c = len(mapped_features) - a
            # not in feature_set doesn't match go_id
            d = len(feature_ids) - len(feature_set_ids) - c

            fisher_value = fisher.pvalue(a, b, c, d)
            if statistical_significance == 'left_tailed':
                raw_p_value = self._round(fisher_value.left_tail)
            elif statistical_significance == 'right_tailed':
                raw_p_value = self._round(fisher_value.right_tail)
            elif statistical_significance == 'two_tailed':
                raw_p_value = self._round(fisher_value.two_tail)
            else:
                raise ValueError('Improper statistical_significance value')

            all_raw_p_value.append(raw_p_value)
            go_info_map.update({go_id: {'raw_p_value': raw_p_value,
                                        'num_in_ref_genome': len(mapped_features),
                                        'num_in_subset_feature_set': a,
                                        'pos': pos,
                                        'mapped_features': mapped_features}})
            pos += 1

        stats = importr('stats')
        adjusted_p_values = stats.p_adjust(FloatVector(all_raw_p_value), method='fdr')

        for go_id, go_info in go_info_map.items():
            if go_id not in ontology_hash:
                continue

            adjusted_p_value = self._round(adjusted_p_values[go_info.get('pos')])
            namespace = ontology_hash[go_id]['namespace']        ############################## prob change below
            enrichment_map.update({go_id: {'raw_p_value': go_info.get('raw_p_value'),
                                           'adjusted_p_value': adjusted_p_value,
                                           'num_in_ref_genome': go_info.get('num_in_ref_genome'),
                                           'num_in_subset_feature_set':
                                           go_info.get('num_in_subset_feature_set'),
                                           'go_term': go_id_go_term_map.get(go_id),
                                           'namespace': namespace.split("_")[1][0].upper(),
                                           'mapped_features': go_info.get('mapped_features')}})

        returnVal = {'result_directory': result_directory}
        report_output = self._generate_report(enrichment_map,
                                              result_directory,
                                              params.get('workspace_name'),
                                              feature_id_go_id_list_map,
                                              feature_set_ids,
                                              genome_ref,
                                              go_id_parent_ids_map,
                                              feature_ids)

        returnVal.update(report_output)


        ###############################################
        # Part of what we want this program to do is map GO terms onto
        # different types of listing methods (EG: EC, KEGG, SEED)
        # Current idea on how to do this is to convert GO terms to EC
        # terms taken from
        # http://current.geneontology.org/ontology/external2go/ec2go
        # and then translate that into SEED data
        # can use EC conversion data:
        # https://github.com/kbaseapps/GenomeFileUtil/blob/e5f3cb9dd1f6fe3d313551d4780c8dca422e9e63/data/cath_ontology_mapping.json
        #or GO conversion data:
        #https://github.com/kbaseapps/GenomeFileUtil/blob/e5f3cb9dd1f6fe3d313551d4780c8dca422e9e63/data/go_ontology_mapping.json
        ###############################################



        return returnVal
