<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Website</title>

    <style>
        body {
            margin: 20px; /* Add some margin to the body to provide space */
            background-color: #F5DEB3; /* Change the background color to light tan*/
        }

        .container {
            max-width: 800px; /* Set maximum width for better readability */
            margin: 0 auto; /* Center the content horizontally */
            padding: 20px; /* Add padding around the content */
        }

        .column {
            margin-bottom: 20px; /* Add margin between columns */
        }

        .greeting {
            font-size: 1.2em; /* Increase the font size of the greeting */
            margin-bottom: 20px; /* Add margin to the bottom of the greeting */
            text-align: justify; /* Justify the text */
        }

        h2 {
            margin-bottom: 10px; /* Add margin below headings */
        }

        ul {
            list-style-type: none; /* Remove bullets from lists */
            padding: 0; /* Remove default padding */
        }

        img {
            display: block;
            width: 100%; /* Make images responsive, occupying full width */
            height: auto;
            margin-bottom: 10px; /* Add margin below images */
        }

        .contactInfo {
            text-align: center; /* Center contact information */
        }
    </style>
</head>
<body>

<div class="container">
    <!-- Second Column (Above First Column) -->
    <div class="column">
        <!-- Image -->
        <img src="https://github.com/ElizabethThompson98/ElizabethThompson98.github.io/blob/main/Directory_Photo.jpg?raw=true" alt="">

        <!-- Contact Information -->
        <div class="contactInfo">
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
    </div>

    <!-- First Column (Below Second Column) -->
    <div class="column">
        <!-- Greeting at the top of the first column -->
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
            <li><a href="https://www.westerncriminology.org/documents/conference_proceedings/WSC_2024_Conference_Program.pdf"> Western Society of Criminology 50th Annual Conference</a> "Perceptions of Force Severity, Operationalizing Reasonableness With a Repeated Measures Design" by David Makin, Elizabeth Thompson, Mary McMillin, Bala Krishnamoorthy, and Dale Willits, Long Beach CA, February 9, 2024</li>
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

        <!-- Papers -->
        <h2>Papers (in progress)</h2>
        <ul>
            <li> <a href="DQE_Final_Draft.pdf" target="_blank">Comparing TDA Methods of Time Series Analysis</a> 
                 <a href="R.Code Main.R" target="_blank">with associated R code</a> </li>
            <li> <a href="Perceptions_of_Force_Update_March_2024.pdf" target="_blank">Predictive Models for Perceptions of Police-Community Interactions Involving Use of Force</a> 
                 <a href="Perceptions_UofF_Code.ipynb" target="_blank"> with associated Python code</a> </li>
            <li> <a href="Mapper_Analysis.pdf" target="_blank">Factors That Impact Uses of Force in Police-Community Interactions</a> 
                 <a href="UofF Mapper.ipynb" target="_blank">with associated Python code </a> </li>
        </ul>
    </div>
</div>

</body>
</html>
