jQuery(document).ready(function() {
  // Global variables to store user answers
  let answer1 = "";
  let answer2 = "";

  // Function to update choices for question 1
  function updateChoices1() {
    // Make an API call to Flask to get choices for question 1
    fetch('/api/question1')
      .then(response => response.json())
      .then(choices => {
        // Clear existing choices
        jQuery('#choices1').empty();

        // Add new choices
        choices.forEach(function(choice) {
          jQuery('#choices1').append('<li><button class="choice">' + choice + '</button></li>');
        });
      })
      .catch(function() {
        alert('Failed to retrieve choices for question 1.');
      });
  }

  // Function to update choices for question 2 based on answer 1
  function updateChoices2() {
    // Make an API call to Flask to get choices for question 2 based on answer 1
    fetch('/api/question2', {
      method: 'POST',
      body: new URLSearchParams({ answer1: answer1 })
    })
      .then(response => response.json())
      .then(choices => {
        // Clear existing choices
        jQuery('#choices2').empty();

        // Add new choices
        choices.forEach(function(choice) {
          jQuery('#choices2').append('<li><button class="choice">' + choice + '</button></li>');
        });
      })
      .catch(function() {
        alert('Failed to retrieve choices for question 2.');
      });
  }

  // Function to update choices for question 3 based on answer 2
  function updateChoices3() {
    // Make an API call to Flask to get choices for question 3 based on answer 2
    fetch('/api/question3', {
      method: 'POST',
      body: new URLSearchParams({ answer2: answer2 })
    })
      .then(response => response.json())
      .then(choices => {
        // Clear existing choices
        jQuery('#choices3').empty();

        // Add new choices
        choices.forEach(function(choice) {
          jQuery('#choices3').append('<li><button class="choice">' + choice + '</button></li>');
        });
      })
      .catch(function() {
        alert('Failed to retrieve choices for question 3.');
      });
  }

  // Event handler for clicking on choices
  jQuery('.choices').on('click', '.choice', function() {
    var choice = jQuery(this).text();

    // Determine which question the choice belongs to
    var questionId = jQuery(this).closest('.question').attr('id');

    // Store user answers and update choices accordingly
    if (questionId === 'question1') {
      answer1 = choice;
      updateChoices2();
    } else if (questionId === 'question2') {
      answer2 = choice;
      updateChoices3();
    } else if (questionId === 'question3') {
      // User has answered all three questions, perform final processing
      // TODO: Make another API call to Flask to handle the final processing
      alert('All questions answered! Processing complete.');
    }
  });

  // Initial setup - Retrieve choices for question 1
  updateChoices1();
});
