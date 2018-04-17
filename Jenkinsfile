pipeline {
  agent {
kubernetes {
      //cloud 'kubernetes'
      label 'mypod'
      containerTemplate {
        name 'dealii'
        image 'ubuntu'
        ttyEnabled true
        command 'cat'
      }
}
  }
  stages {
    stage('astyle') {
      steps {
        container('dealii'){
        sh 'id'
        sh 'ls -la'
        sh './doc/indent'
        sh 'git diff | tee astyle-changes.diff'
        archiveArtifacts artifacts: 'astyle-changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
        }
      }
    }
    stage('build-gcc-fast') {
      steps {
        container('dealii'){
        sh '''
          mkdir -p build-gcc-fast
          cd build-gcc-fast
          cmake -G "Ninja" gcc -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_USE_PETSC=OFF -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON ..
          ninja
        '''
        }
      } 
    }
    stage('test') {
      steps {
        container('dealii'){
        sh 'cd build-gcc-fast && ../cmake/generate_reference_output.sh'
        archiveArtifacts artifacts: 'changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
        }
      }
    }
  }
}
