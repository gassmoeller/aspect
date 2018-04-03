pipeline {
  agent {
      docker { image 'gassmoeller/aspect-tester:astyle' }
  }
  stages {
    stage('astyle') {
      steps {
        sh './doc/indent'
        sh 'git diff | tee logfile'
        archiveArtifacts artifacts: 'logfile', fingerprint: true
        sh 'git diff --exit-code --name-only'
      }
    }
  }
}
