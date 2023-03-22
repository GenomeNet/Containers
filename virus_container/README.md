Container to perform prediction of viruses (binary / genus level supported)

## About

- Input: FASTA files
- Output: file with predictions written to the /output folder of the mounted folder

## Create

```
sudo docker build -t genomenet/virus .
sudo docker push genomenet/virus
podman pull genomenet/virus
```

Then use the python scrip tin the root folder to run the container. 