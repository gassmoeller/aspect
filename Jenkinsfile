#!groovy

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
    stage("Debug") {
      steps {
          echo "Running build ${env.BUILD_ID} on ${env.NODE_NAME}, env=${env.NODE_ENV}"
          sh 'printenv'
          sh 'echo $CHANGE_AUTHOR_EMAIL'
          sh "echo ${env.CHANGE_AUTHOR_EMAIL}"
      }
    }

    stage ("Check execution") {
      when {
	expression { return ! ( "${env.TRUST_BUILD}" == "true" || "${env.CHANGE_BRANCH}" == "jenkins-pipeline" || ["heister@clemson.edu", "rene.gassmoeller@mailbox.org"].contains("${env.CHANGE_AUTHOR_EMAIL}") )
	    }
      }
      steps {
	  echo "Please ask an admin to rerun jenkins with TRUST_BUILD=true"
	    sh "exit 1"
      }
    }

    stage('Check indentation') {
      steps {
        sh './doc/indent'
        sh 'git diff | tee astyle-changes.diff'
        archiveArtifacts artifacts: 'astyle-changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
        }
    }
    stage('Build') {
      steps {
        sh '''
          mkdir -p build-gcc-fast
          cd build-gcc-fast
          cmake -G "Ninja" gcc -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_USE_PETSC=OFF -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON ..
          ninja
        '''
      } 
    }
    stage('Run tests') {
      steps {
        sh '''
          cd build-gcc-fast
          ASPECT_TESTS_VERBOSE=1 ../cmake/generate_reference_output.sh
        '''
        archiveArtifacts artifacts: 'build-gcc-fast/changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
      }
    }
  }
  post {
    always {
      deleteDir()
    }
  }
}
