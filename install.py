#!/usr/bin/env python3

import subprocess
import platform

# Install necessary libraries
libraries = ["argparse", "pathlib2", "termcolor"]
subprocess.run(["pip3", "install"] + libraries, check=True)
print(f"Installed libraries: {libraries}")

# Check if Podman is installed, and install it if it is not present
if subprocess.run(["which", "podman"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0:
    # Prompt the user to confirm installing Podman
    user_input = input("Podman is not installed. Do you want to install it? (y/n) ").lower()
    while user_input not in ["y", "n"]:
        user_input = input("Please enter 'y' to install Podman or 'n' to exit: ").lower()
    if user_input == "n":
        print("Aborting installation.")
        exit()

    # Install Podman based on the current distribution
    print("Installing Podman...")
    distribution = platform.linux_distribution()[0]
    if distribution == "Ubuntu":
        subprocess.run(["sudo", "apt", "update"], check=True)
        subprocess.run(["sudo", "apt", "install", "-y", "podman"], check=True)
    elif distribution == "CentOS Linux":
        subprocess.run(["sudo", "yum", "update"], check=True)
        subprocess.run(["sudo", "yum", "install", "-y", "podman"], check=True)
    print("Podman installed successfully.")
else:
    print("Podman is already installed.")
