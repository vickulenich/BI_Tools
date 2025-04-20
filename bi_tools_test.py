import pytest
import os
from Bio import SeqIO
from bi_tools import filter_fastq

@pytest.fixture
def create_temp_fastq_file(tmp_path):
    """Create a temporary FASTQ file for testing."""
    fastq_content = """@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%+++**))**))**))**))**))**))**))**))**))**))"""
    fastq_file = tmp_path/"test.fastq"
    with open(fastq_file, "w") as f:
        f.write(fastq_content)
    return fastq_file

def test_filter_fastq_valid(create_temp_fastq_file):
    input_file = create_temp_fastq_file()
    output_file = create_temp_fastq_file.parent/"output.fastq"
    filter_fastq(str(input_file), str(output_file), (0, 100), (0, 30), 20)
    
    assert os.path.exists(output_file)
    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 1

def test_filter_fastq_gc_bounds(create_temp_fastq_file):
    input_file = create_temp_fastq_file
    output_file = create_temp_fastq_file.parent/"output_gc.fastq"
    filter_fastq(str(input_file), str(output_file), (50, 60), (0, 30), 20)
    
    assert os.path.exists(output_file)
    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 0

def test_filter_fastq_length_bounds(create_temp_fastq_file):
    input_file = create_temp_fastq_file
    output_file = create_temp_fastq_file.parent/"output_length.fastq"
    filter_fastq(str(input_file), str(output_file), (0, 100), (10, 20), 20)
    
    assert os.path.exists(output_file)
    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 0

def test_filter_fastq_quality_threshold(create_temp_fastq_file):
    input_file = create_temp_fastq_file
    output_file = create_temp_fastq_file.parent/"output_quality.fastq"
    filter_fastq(str(input_file), str(output_file), (0, 100), (0, 30), 30)
    
    assert os.path.exists(output_file)
    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 0

def test_filter_fastq_file_not_found():
    output_file = "output.fastq"
    result = filter_fastq("non_existent.fastq", output_file)
    
    assert result is None

def test_filter_fastq_empty_file(tmp_path):
    """Test filtering on an empty FASTQ file"""
    input_file = tmp_path / "empty.fastq"
    output_file = tmp_path / "output_empty.fastq"
    input_file.touch()
    filter_fastq(str(input_file), str(output_file), (0, 100), (0, 30), 20)
    
    assert os.path.exists(output_file)
    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 0

def test_filter_fastq_invalid_quality_scores(create_temp_fastq_file):
    input_file = create_temp_fastq_file
    output_file = create_temp_fastq_file.parent/"output_invalid_quality.fastq"
    filter_fastq(str(input_file), str(output_file), (0, 100), (0, 30), -1)
    
    assert os.path.exists(output_file)
    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 1

def test_filter_fastq_gc_content_edge_cases(create_temp_fastq_file):
    input_file = create_temp_fastq_file
    output_file = create_temp_fastq_file.parent/"output_gc_edge.fastq"
    filter_fastq(str(input_file), str(output_file), (0, 100), (0, 30), 20)
    
    assert os.path.exists(output_file)
    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 1

def test_filter_fastq_length_edge_cases(create_temp_fastq_file):
    input_file = create_temp_fastq_file
    output_file = create_temp_fastq_file.parent/"output_length_edge.fastq"
    filter_fastq(str(input_file), str(output_file), (0, 100), (15, 30), 20)
    
    assert os.path.exists(output_file)
    records = list(SeqIO.parse(str(output_file), "fastq"))
    assert len(records) == 1