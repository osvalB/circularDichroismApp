#!/bin/bash

printMessageHowToTestDocker () {

    echo -e "\nTo test $1 in the developer machine run $ sudo docker run -p 3838:3838 $2"
    echo -e  "and then open http://localhost:3838/ in your browser"

}

echo -e "\nBuilding Image for Circular dichroism..."
docker build -t chirakit_rafael -f ./Dockerfile_circularDichroism .

printMessageHowToTestDocker 'ChiraKit' 'chirakit_rafael'
