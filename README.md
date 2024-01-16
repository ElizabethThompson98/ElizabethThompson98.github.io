<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ElizabethThompson98</title>

    <style>
        body {
            margin: 20px; /* Add some margin to the body to provide space */
            display: flex;
            flex-wrap: wrap; /* Allow items to wrap to the next line if needed */
        }

        .column {
            flex: 1; /* Each column takes up equal space */
            max-width: 50%; /* Limit each column to half the width of the page */
            margin-right: 0; /* Remove the margin between columns */
            box-sizing: border-box; /* Include padding and border in the total width */
        }

        .greeting {
            font-size: 1.2em; /* Increase the font size of the greeting */
            margin-bottom: 20px; /* Add margin to the bottom of the greeting */
        }      

        p {
            text-align: justify; /* Justify the text */
        }

        h2 {
            clear: right; /* Clear the float to ensure headings start below the image */
        }

        /* Style to remove bullets from all lists */
        ul {
            list-style-type: none;
            padding: 0; /* Remove default padding */
        }
    </style>
</head>
<body>

    <div class="column">
        <!-- Greeting at the top of the left column -->
        <p class="greeting">
            Hi there! My name is Elizabeth Thompson. I am a research and teaching assistant at Washington State University. I love teaching Calculus, and my research is in Topological Data Analysis and Machine Learning!
            <!-- Link to the PDF file -->
            <a href="Graduate CV.pdf" target="_blank">Here is a link to my CV.</a>
        </p>

        <!-- Education -->
        <h2>Education</h2>
        <ul>
            <li>PhD in Mathematics, Washington State University (WSU), 2021-Present</li>
            <li>BS in Mathematics & Secondary Education, Linfield University, 2017-2021</li>
        </ul>

        <!-- Presentations -->
        <h2>Presentations</h2>
        <ul>
            <li><a href="https://www.youtube.com/watch?v=iu7vNOukW00&t=3550s">What Can Donuts Tell Us About Body Cam Data?</a> 2nd Place Three Minute Thesis Finals, WSU, March 29, 2023</li>        
            <li><a href="https://georgefoxuniversity.regfox.com/pacific-northwest-section-of-the-mathematical-association-of-america"> PNW MAA Meeting</a> "Motivating Student Learning with Mathematical Modeling in Calculus for Life Science" by Alexander Dimitrov and Elizabeth Thompson, George Fox University, March 18, 2023</li>
        </ul>

        <!-- Teaching -->
        <h2>Teaching</h2>
        <ul>
            <li>Math 202: Calculus for Business & Economics Lecture, Fall 2023</li>
            <li>Math 103: Algebra Methods & Functions Lecture, Summer 2023</li>
            <li>Math 108: Trigonometry Lecture, Summer 2022 & Summer 2023</li>
            <li>Math 140 Lab: Calculus for Life Scientists, Fall 2023-Spring 2023</li>
            <li>Math 171 Global Campus: Calculus I Lecture & Lab, Summer 2022</li>
        </ul>
    </div>

    <div class="column">
        <!-- Image and Contact Information -->
        <img src="https://github.com/ElizabethThompson98/ElizabethThompson98.github.io/blob/main/Directory_Photo.jpg?raw=true" alt="Elizabeth Thompson's Photo" style="max-width: 100%; height: auto;">

        <p>
            <strong>Contact Information</strong><br>
            Email: elizabeth.thompson1@wsu.edu<br>
            Office Address:<br>
            Vancouver Undergraduate Building<br>
            Room 251<br>
            Washington State University<br>
            Vancouver WA, 98661
        </p>
    </div>

</body>
</html>
