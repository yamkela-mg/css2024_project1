#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 23:11:28 2024

@author: yamkelamgwatyu
"""

import pandas as pd

df = pd.read_csv("/Users/yamkelamgwatyu/Documents/CSS2024_training/movie_dataset.csv")

print(df.info())

#Q1
highest_rated_movie = df.loc[df['Rating'].idxmax()]
print(highest_rated_movie)

#Q2
average_revenue = df['Revenue (Millions)'].mean()
print("Average Revenue of All Movies:", average_revenue)

#Q3
filtered_df = df[(df['Year'] >= 2015) & (df['Year'] <= 2017)]
average_revenue_2015_2017 = filtered_df['Revenue (Millions)'].mean()
print("Average Revenue of Movies from 2015 to 2017:", average_revenue_2015_2017)

#Q4
movies_2016 = df[df['Year'] == 2016]
count_movies_2016 = len(movies_2016)
print("Number of Movies Released in 2016:", count_movies_2016)

#Q5
nolan_movies = df[df['Director'] == 'Christopher Nolan']
count_nolan_movies = len(nolan_movies)
print("Number of Movies Directed by Christopher Nolan:", count_nolan_movies)

#Q6
highly_rated_movies = df[df['Rating'] >= 8.0]
count_highly_rated_movies = len(highly_rated_movies)
count_highly_rated_movies = len(highly_rated_movies)
print("Number of Movies with a Rating of at Least 8.0:", count_highly_rated_movies)

#Q7
nolan_movies = df[df['Director'] == 'Christopher Nolan']
median_rating_nolan_movies = nolan_movies['Rating'].median()
print("Median Rating of Movies Directed by Christopher Nolan:", median_rating_nolan_movies)

#Q8
average_rating_by_year = df.groupby('Year')['Rating'].mean()
year_highest_avg_rating = average_rating_by_year.idxmax()
highest_avg_rating = average_rating_by_year.max()
print("Year with the Highest Average Rating:", year_highest_avg_rating)
print("Highest Average Rating:", highest_avg_rating)

#Q9
movies_2006 = df[df['Year'] == 2006]
movies_2016 = df[df['Year'] == 2016]
num_movies_2006 = len(movies_2006)
num_movies_2016 = len(movies_2016)
percentage_increase = ((num_movies_2016 - num_movies_2006) / num_movies_2006) * 100
print("Percentage Increase in Number of Movies (2006 to 2016):", percentage_increase)

#Q10
all_actors = df['Actors'].str.split(',').explode().str.strip()
most_common_actor = all_actors.mode().iloc[0]
print("Most Common Actor in All Movies:", most_common_actor)

#Q11
all_genres = df['Genre'].str.split(',').explode().str.strip()
num_unique_genres = all_genres.nunique()
print("Number of Unique Genres:", num_unique_genres)

#Q12
numerical_columns = ['Year', 'Runtime (Minutes)', 'Rating', 'Votes', 'Revenue (Millions)']
numerical_df = df[numerical_columns]
correlation_matrix = numerical_df.corr()
print("Correlation Matrix:")
print(correlation_matrix)

