#!/bin/bash
echo "Stopping vol2bird container ..."
docker kill vol2bird
echo "Removing stopped vol2bird container ..."
docker rm vol2bird
