pipeline {
  agent {
      docker { image 'gassmoeller/aspect-tester:astyle' }
  }
  stages {
    stage('astyle') {
      steps {
        sh './doc/indent'
        sh 'git diff | tee astyle-changes.diff'
        archiveArtifacts artifacts: 'astyle-changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
      }
    }
    stage('build') {
      steps {
        sh 'cmake -G "Ninja" gcc -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_USE_PETSC=off -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON .'
        sh 'ninja'
      } 
    }
    stage('test') {
      steps {
        sh 'ctest --output-on-failure -j4 || { echo "test FAILED"; }'
        sh 'ninja generate_reference_output'
      }
    }
  }
}
