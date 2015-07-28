__author__ = 'rileyd'


def write_to_xls(workbook, molecules):
    """"
    Given an array of JSON  objects, write them out to an upload ready XLSX
    File. Column headers will have the full attribute path of the data,
    on record will be stored per line.
    **mandatory child objects will be denormalized into the same row as their
    parent**
    """
    feat_sheet = workbook.add_worksheet('dnafeatures')
    seen_feat = {}
    ft_row = 0
    ft_cols = ['category.name', 'accession', 'name', 'pattern.bases',
               'pattern.sha1', 'description', 'properties']
    col = 0
    # Header
    for col, heading in enumerate(ft_cols):
        feat_sheet.write_string(ft_row, col, heading)
        col += 1
        feat_sheet.write_string(ft_row, col, "Extracted From")

    ft_row += 1

    for mol in molecules:
        annotations = mol.get('dnafeatures')
        for ann in annotations:
            ft = ann.get('dnafeature')
            ft_accession = ft.get('accession')
            definition = "{mol_accession}:{start}..{end}".format(
                mol_accession=mol['accession'],
                start=ann['start'],
                end=ann['end']
            )
            if ft_accession in seen_feat:
                # now check the sequence
                seen_bases = seen_feat[ft_accession]['pattern']['sha1']
                new_bases = ft['pattern']['sha1']
                if new_bases == seen_bases:
                    log.debug("Skipping already-seen {0}".format(ft_accession))
                    continue
                else:
                    log.warn("{0} re-defined with new nucleotide pattern in {1}".format(
                        ft_accession, definition))

            for col, attr in enumerate(ft_cols):
                parent = ft
                if '.' in attr:
                    attr_path = attr.split('.')
                    for child_name in attr_path:
                        parent = parent[child_name]
                    value = parent
                else:
                    value = ft.get(attr)
                if isinstance(value, basestring):
                    feat_sheet.write_string(ft_row, col, value)
                elif isinstance(value, dict):
                    if value == {}:
                        continue
                    flat_value = json.dumps(value)
                    feat_sheet.write_string(ft_row, col, flat_value)
                else:
                    feat_sheet.write(ft_row, col, value)

            # Extraction Source
            feat_sheet.write_string(ft_row, (col + 1), definition)
            ft_row += 1
            seen_feat[ft_accession] = ft  # stash the feature

