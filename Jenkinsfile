pipeline {
  agent {
    docker {
        image 'dealii/dealii:v8.5.1-gcc-mpi-fulldepscandi-debugrelease'
	label 'has-docker'
	args '-v /home/docker/jenkins:/home/dealii/jenkins'
    }
  }

  parameters {
    booleanParam(defaultValue: false, description: 'do we trust this user to run the testsuite?', name: 'TRUST_BUILD')
  }

  stages {
    stage("pre") {
      steps {
          echo "Running build ${env.BUILD_ID} on ${env.NODE_NAME}, env=${env.NODE_ENV}"
          sh 'printenv'
          echo '${params.TRUST_BUILD}'
          echo '${TRUST_BUILD}'
          echo '${env.TRUST_BUILD}'
      }
    }

    stage('astyle') {
      steps {
        sh './doc/indent'
        sh 'git diff | tee astyle-changes.diff'
        archiveArtifacts artifacts: 'astyle-changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
        }
    }
    stage('build-gcc-fast') {
      steps {
        sh '''
          mkdir -p build-gcc-fast
          cd build-gcc-fast
          cmake -G "Ninja" gcc -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_USE_PETSC=OFF -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON ..
          ninja
        '''
      } 
    }
    stage('test') {
      steps {
        sh '''
          sed -i 's/mpirun/mpirun --allow-run-as-root/' tests/CMakeLists.txt
          git add tests/CMakeLists.txt && git commit -m 'tester specific changes' --author 'tester <tester@tester.com>'
          cd build-gcc-fast
          ASPECT_TESTS_VERBOSE=1 ../cmake/generate_reference_output.sh
        '''
        archiveArtifacts artifacts: 'build-gcc-fast/changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
      }
    }
  }
}
