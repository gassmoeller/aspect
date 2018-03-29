pipeline {
  agent {
      docker { image 'gassmoeller/aspect-tester:astyle' }
  }
  stages {
    stage('astyle') {
      steps {
        echo "I am alive"
        sh 'pwd'
        sh 'whoami'
        sh 'ls /home/dealii/script.sh'
      }
    }
  }
}
