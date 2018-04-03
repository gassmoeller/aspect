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
        sh 'ls'
        sh 'cmake .'
        sh 'make -j'
      } 
    }
  }
}
