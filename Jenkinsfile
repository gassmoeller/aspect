pipeline {
  agent {
kubernetes {
      //cloud 'kubernetes'
      label 'mypod'
      containerTemplate {
        name 'aspect-tester'
        image 'gassmoeller/aspect-tester:stable2'
        ttyEnabled true
        command 'cat'
      }
}
  }
  stages {
    stage('astyle') {
      steps {
        container('aspect-tester'){
        sh 'ls -la'
        sh 'ls -la /home/jenkins/workspace/'
        sh 'pwd'
        sh 'id'
        sh 'su - jenkins && id && exit'
        sh 'cd && mkdir aspect-build'
        sh './doc/indent'
        sh 'git diff | tee astyle-changes.diff'
        archiveArtifacts artifacts: 'astyle-changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
        }
      }
    }
    stage('build-gcc-fast') {
      steps {
        container('aspect-tester'){
        sh '''
          su jenkins
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
        container('aspect-tester'){
        sh 'su jenkins && cd build-gcc-fast && ASPECT_TESTS_VERBOSE=1 ../cmake/generate_reference_output.sh'
        archiveArtifacts artifacts: 'build-gcc-fast/changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
        }
      }
    }
  }
}
