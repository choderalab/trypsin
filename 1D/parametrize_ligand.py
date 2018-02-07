# coding=utf-8
"""
This script parametrizes 1D
"""

from sys import argv
from protons.app.ligands import generate_protons_ffxml, gaff_default
from protons.app import logger
from protons.app.logger import log
from lxml import etree


log.setLevel(logger.logging.DEBUG)


def parse_args(args):

    parsed_args = None
    # no arguments
    if len(args) == 1:
        parsed_args= dict(infile =  "1D_allH.mae", outfile="1D.xml", hydrogens="1D-hydrogens.xml")
    # 1 argument, presumed to be the input file
    elif len(args) == 2:
        # name output similar to input
        name = args[1].split('.')[0]
        parsed_args = dict(infile=args[1], outfile=name+".xml", hydrogens=name+"-hydrogens.xml")
    # 2 arguments, first is input second is output
    elif len(args) == 3:
        name = args[2].split('.')[0]
        parsed_args= dict(infile=args[1], outfile=args[2], hydrogens=name+"-hydrogens.xml")
    # wrong cases
    else:
        raise ValueError("Please provide input file as the first argument, and optionally an output name as the second argument")

    return parsed_args


def create_hydrogen_definitions(inputfile, outputfile, gaff=gaff_default):
    """
    Generates hydrogen definitions for a small molecule residue template.
    
    Parameters
    ----------
    inputfile - str
        a forcefield XML file defined using Gaff atom types
    outputfile - str
        Name for the XML output file
    """

    gafftree = etree.parse(gaff, etree.XMLParser(remove_blank_text=True, remove_comments=True))
    xmltree = etree.parse(inputfile, etree.XMLParser(remove_blank_text=True, remove_comments=True))
    # Output tree
    hydrogen_definitions_tree = etree.fromstring('<Residues/>')
    hydrogen_types = _find_hydrogen_types(gafftree)

    for residue in xmltree.xpath('Residues/Residue'):
        hydrogen_file_residue = etree.fromstring("<Residue/>")
        hydrogen_file_residue.set('name', residue.get('name'))
        # enumerate hydrogens in this list
        hydrogens = list()
        # Loop through atoms to find all hydrogens
        for atom in residue.xpath('Atom'):
            if atom.get('type') in hydrogen_types:
                # Find the parent atom
                for bond in residue.xpath('Bond'):
                    atomname1 = bond.get('atomName1')
                    atomname2 = bond.get('atomName2')
                    # There should be only one bond containing this hydrogen
                    if atom.get('name') == atomname1:
                        hydrogens.append(tuple([atomname1, atomname2]))
                        break
                    elif atom.get('name') == atomname2:
                        hydrogens.append(tuple([atomname2, atomname1]))
                        break

        # Loop through all hydrogens, and create definitions
        for name, parent in hydrogens:
            h_xml = etree.fromstring("<H/>")
            h_xml.set("name", name)
            h_xml.set("parent", parent)
            hydrogen_file_residue.append(h_xml)
        hydrogen_definitions_tree.append(hydrogen_file_residue)
    # Write output
    xmlstring = etree.tostring(hydrogen_definitions_tree, encoding="utf-8", pretty_print=True, xml_declaration=False)
    xmlstring = xmlstring.decode("utf-8")
    with open(outputfile, 'w') as fstream:
        fstream.write(xmlstring)


def _find_hydrogen_types(gafftree):
    """
    Find all atom types that describe hydrogen atoms.
    
    Parameters
    ----------
    gafftree - lxml.ElementTree
        A GAFF input xml file that contains atom type definitions.

    Returns
    -------
    set - names of all atom types that correspond to hydrogen

    """

    # Detect all hydrogen types by element and store them in a set
    hydrogen_types = set()
    for atomtype in gafftree.xpath('AtomTypes/Type'):
        if atomtype.get('element') == "H":
            hydrogen_types.add(atomtype.get('name'))

    return hydrogen_types


def main(args):
    """
    Run the program
    Parameters
    ----------
    args - cmd line arguments

    """
    import os
    parsed_args = parse_args(args)
    epikoutputdir='epik_output'
    if not os.path.isdir(epikoutputdir):
        os.makedirs(epikoutputdir)

    result = generate_protons_ffxml(parsed_args['infile'], parsed_args['outfile'], tmpdir=epikoutputdir, remove_temp_files=False, pH=7.4, resname="STI")
    create_hydrogen_definitions(parsed_args['outfile'], parsed_args['hydrogens'])

if __name__ == "__main__":
    main(argv)
