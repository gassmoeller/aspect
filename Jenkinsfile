pipeline {
  agent {
      docker { image 'dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease' }
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
        sh 'cd build-gcc-fast && ../cmake/generate_reference_output.sh'
        archiveArtifacts artifacts: 'changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
      }
    }
  }
  post {
        always {
            deleteDir() /* clean up our workspace */
        }
  }
}
