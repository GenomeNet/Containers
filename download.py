#!/usr/bin/env python3

import argparse
import subprocess

def is_podman_installed():
    try:
        subprocess.run(["podman", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except FileNotFoundError:
        return False

def download_container(container_name):
    if not is_podman_installed():
        print("Error: Podman is not installed.")
        print("Please follow the installation instructions at https://podman.io/getting-started/installation.html")
        return

    try:
        subprocess.run(["podman", "pull", f"docker.io/{container_name}"], check=True)
        print(f"Container '{container_name}' downloaded successfully.")
    except subprocess.CalledProcessError:
        print(f"Error: Failed to download container '{container_name}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download a container using podman from docker.io")
    parser.add_argument("--container", required=True, help="Name of the container to download from docker.io")

    args = parser.parse_args()
    download_container(args.container)
