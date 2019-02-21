#!/bin/sh

set -ev

# Run all debug tests
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  echo "Hello from Travis"
  pytest
fi
