import requests
import time
import os
from queue import Queue
from threading import Thread

MAX_CONCURRENT_JOBS = 3  # Limit on simultaneous submissions
MAX_RETRIES = 3  # Maximum retries for failed jobs

def split_fasta(input_file, batch_size):
    """Split a FASTA file into smaller files with a specified number of sequences."""
    with open(input_file, "r") as infile:
        batches = []
        batch = []
        for line in infile:
            if line.startswith(">") and len(batch) >= batch_size:
                batches.append(batch)
                batch = []
            batch.append(line)
        if batch:
            batches.append(batch)
    
    # Write batches to separate files
    batch_files = []
    for i, batch in enumerate(batches):
        batch_file = f"batch_{i + 1}.fasta"
        with open(batch_file, "w") as outfile:
            outfile.writelines(batch)
        batch_files.append(batch_file)
    return batch_files

def submit_to_interpro(batch_file):
    """Submit a batch of sequences to InterProScan and return the job ID."""
    url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run/"
    with open(batch_file, "r") as fasta:
        data = {
            "email": "your_email@example.com",  # Replace with your email
            "sequence": fasta.read()
        }
        response = requests.post(url, data=data)
        if response.status_code == 200:
            job_id = response.text.strip()
            print(f"Submitted batch {batch_file}. Job ID: {job_id}")
            return job_id
        else:
            print(f"Failed to submit batch {batch_file}: {response.text}")
            return None

def check_status(job_id):
    """Check the status of a submitted job."""
    url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
    response = requests.get(url)
    return response.text.strip()

def fetch_results(job_id, output_file):
    """Fetch the results of a completed job and save them to a TSV file."""
    url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/tsv"
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_file, "w") as outfile:
            outfile.write(response.text)
        print(f"Results saved to {output_file}")
        return True
    else:
        print(f"Failed to fetch results for Job ID {job_id}: {response.text}")
        return False

def process_batches(batch_files, final_output):
    """Submit and process batches with retry logic and concurrent job handling."""
    job_queue = Queue()
    completed_jobs = []
    retries = {}

    def worker():
        while not job_queue.empty():
            batch_file = job_queue.get()
            retry_count = retries.get(batch_file, 0)
            job_id = submit_to_interpro(batch_file)
            if job_id:
                while True:
                    status = check_status(job_id)
                    print(f"Job {job_id} status: {status}")
                    if status == "FINISHED":
                        output_file = f"{job_id}.tsv"
                        if fetch_results(job_id, output_file):
                            completed_jobs.append(output_file)
                        break
                    elif status == "FAILURE":
                        print(f"Job {job_id} failed!")
                        if retry_count < MAX_RETRIES:
                            print(f"Retrying batch {batch_file} ({retry_count + 1}/{MAX_RETRIES})")
                            retries[batch_file] = retry_count + 1
                            job_queue.put(batch_file)
                        break
                    time.sleep(30)  # Check status every 30 seconds
            else:
                print(f"Submission failed for batch {batch_file}")
            job_queue.task_done()

    # Add batches to the queue
    for batch_file in batch_files:
        job_queue.put(batch_file)

    # Start worker threads
    threads = []
    for _ in range(min(MAX_CONCURRENT_JOBS, len(batch_files))):
        thread = Thread(target=worker)
        thread.start()
        threads.append(thread)

    # Wait for all jobs to finish
    for thread in threads:
        thread.join()

    # Merge results
    merge_tsv(completed_jobs, final_output)

def merge_tsv(output_files, final_output):
    """Merge multiple TSV files into a single file."""
    with open(final_output, "w") as outfile:
        for i, file in enumerate(output_files):
            with open(file, "r") as infile:
                if i == 0:
                    outfile.write(infile.read())
                else:
                    next(infile)  # Skip the header line
                    outfile.write(infile.read())
    print(f"Merged results into {final_output}")

def main(input_fasta, batch_size, final_output):
    # Split the input FASTA file into smaller batches
    batch_files = split_fasta(input_fasta, batch_size)
    # Process the batches
    process_batches(batch_files, final_output)
    # Clean up temporary files
    for batch in batch_files:
        os.remove(batch)

if __name__ == "__main__":
    input_fasta = "/home/luke/Documents/data/E.faecium_C68/genomic-features/PGAP/annot.faa"  # Replace with your input FASTA file
    batch_size = 99  # Number of sequences per batch
    final_output = "/home/luke/Documents/data/E.faecium_C68/MTase/MTase-search/intrpro_results.tsv"  # Name of the final merged TSV file
    main(input_fasta, batch_size, final_output)
