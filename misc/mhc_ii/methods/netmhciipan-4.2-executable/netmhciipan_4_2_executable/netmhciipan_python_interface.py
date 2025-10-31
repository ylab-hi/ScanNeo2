import os
import pkg_resources
from collections import namedtuple, defaultdict
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

from allele_info import is_user_defined_allele
from iedbtools_utilities.sequence_io import SequenceOutput

import logging
import shutil

logger = logging.getLogger(__name__)
PERL_PATH = shutil.which('perl') or 'perl'

_executable_filename = \
    pkg_resources.resource_filename('netmhciipan_4_2_executable', 'netMHCIIpan')
EXECUTABLE_PATH = os.path.join(os.path.dirname(__file__), _executable_filename)

PredictionInput = namedtuple('PredictionInput', ['sequence', 'allele', 'binding_length'])


def single_prediction(sequence_list, allele_length_2tuple_list, delete_tempfiles=False, el=False):
    """ | *brief*: Provides python interface to netmhciipan 3rd-party tool.
        Performs a prediction for every sequence / allele combination.
        Returns a dictionary { (sequence, allele): [score, score, ...], ... }
        with one entry per prediction made.
    """
    sequence_filepath = sequence_list_to_fasta_file(sequence_list)
    to_delete = [sequence_filepath]

    try:
        all_scores = {}
        # for allele_name in allele_list:
        for allele_name, binding_length in allele_length_2tuple_list:
            if is_user_defined_allele(allele_name):
                # A user-defined allele must be passed via a fasta file.
                t = NamedTemporaryFile(delete=False, mode='w')
                fasta = SequenceOutput.to_fasta(allele_name)
                t.write(fasta)
                t.close()

                user_defined_allele_filepath = sequence_list_to_fasta_file([allele_name])
                cmd = [EXECUTABLE_PATH, '-f', sequence_filepath, '-hlaseq', user_defined_allele_filepath, '-length', str(binding_length)]
                if not el:
                    cmd += ['-BA']
            else:
                # An allele name can be passed as a command-line argument.
                # netmhciipan uses allele names with the '*', ':' and/or '_' stripped
                stripped_allele_name = allele_name.replace('*', '_').replace(':', '')
                if stripped_allele_name.startswith("H2") or stripped_allele_name.startswith("H-2"):
                    stripped_allele_name = stripped_allele_name.replace("H2", "H-2")
                elif not allele_name.startswith("DRB"):
                    stripped_allele_name = "HLA-%s" % stripped_allele_name.replace('_', '')
                cmd = [EXECUTABLE_PATH, '-f', sequence_filepath, '-a', stripped_allele_name, '-length', str(binding_length)]
                if not el:
                    cmd += ['-BA']

            logger.info('Executing: "{}"'.format(' '.join(cmd)))

            get_scores_from_output = get_allele_name_scores_from_netmhciipan_output
            # process = Popen(cmd, stdout=PIPE)
            process = Popen([PERL_PATH] + cmd, stdout=PIPE)
            stdoutdata, stderrdata_ignored = process.communicate()
            stdoutdata = stdoutdata.decode()
            logger.debug('Raw output:\n{}'.format(stdoutdata))

            scores_by_sequence_idx = get_scores_from_output(stdoutdata, el=el)
            # Wrap the scores up nicely to return to the caller
            for seqidx, scores in scores_by_sequence_idx.items():
                sequence = sequence_list[seqidx]
                prediction_input = PredictionInput(sequence, allele_name, binding_length)
                all_scores[prediction_input] = scores
    finally:
        if delete_tempfiles:
            for filepath in to_delete:
                os.unlink(filepath)
    return all_scores


def sequence_list_to_fasta_file(sequence_list):
    """ Writes the sequences in *sequence_list* as fasta sequences to a file
        and returns the filepath.
    """
    t = NamedTemporaryFile(delete=False, mode='w')
    for i, sequence in enumerate(sequence_list):
        t.write('>seq-{}\n'.format(i))
        t.write(sequence)
        t.write('\n')
    t.close()
    return t.name


def get_allele_name_scores_from_netmhciipan_output(netmhciipan_output, el=False):
    """ | *brief*: Parses the string *netmhciipan_output* for scores from an allele name
        |    prediction request.
    """
    logger.debug('netmhciipan_output: %s' % netmhciipan_output )
    netmhciipan_output = [line.split() for line in netmhciipan_output.split('\n') if not line.startswith('#') and line.split() and line.split()[0].isdigit()]
    scores = defaultdict(list)
    for row in netmhciipan_output:
        logger.debug('row: %s' % row)
        # the one which is Identity (for exampole, 'seq-0')
        idx = 6
        sequence_idx = int(row[idx].split('-')[1])
        core = row[4]
        # EL
        if el:
            ic50_score = float(row[idx+1])
        # BA
        else: 
            ic50_score = float(row[idx+5])
        core_score_tuple = (core, ic50_score)
        scores[sequence_idx].append(core_score_tuple)
    return scores

# Retrieved using option -list on the executable
# process = Popen([EXECUTABLE_PATH, '-list'], stdout=PIPE)
process = Popen([PERL_PATH, EXECUTABLE_PATH, '-list'], stdout=PIPE)
allowed_allele_names = process.communicate()[0]

