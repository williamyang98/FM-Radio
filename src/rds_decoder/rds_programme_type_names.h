#pragma once

struct rds_programme_type_name_t {
    const char* description;    // Document description
    const char* display_short;  // 8 characters
    const char* display_long;   // 16 characters
};

// ANNEX F: Programme type codes
// Table F.1: Programme type codes
constexpr int TOTAL_PROGRAMME_TYPES = 32;
const rds_programme_type_name_t PROGRAMME_TYPES[TOTAL_PROGRAMME_TYPES] = {
    { "No programme type or undefined", "None", "None"},
    { "News", "News", "News" },
    { "Current Affairs", "Affairs", "Current Affairs" },
    { "Information", "Info", "Information" },
    { "Sport", "Sport", "Sport" },
    { "Education", "Educate", "Education" },
    { "Drama", "Drama", "Drama" },
    { "Culture", "Culture", "Cultures" },
    { "Science", "Science", "Science" },
    { "Varied", "Varied", "Varied Speech" },
    { "Pop Music", "Pop M", "Pop Music" },
    { "Rock Music", "Rock M", "Rock Music" },
    { "Easy Listening Music", "Easy M", "Easy Listening" },
    { "Light classical", "Light M", "Light Classics M" },
    { "Serious classical", "Classics", "Serious Classics" },
    { "Other Music", "Other M", "Other Music" },
    { "Weather", "Weather", "Weather & Metr" },
    { "Finance", "Finance", "Finance" },
    { "Children’s programmes", "Children", "Children’s Progs" },
    { "Social Affairs", "Social", "Social Affairs" },
    { "Religion", "Religion", "Religion" },
    { "Phone In", "Phone In", "Phone In" },
    { "Travel", "Travel", "Travel & Touring" },
    { "Leisure", "Leisure", "Leisure & Hobby" },
    { "Jazz Music", "Jazz", "Jazz Music" },
    { "Country Music", "Country", "Country Music" },
    { "National Music", "Nation M", "National Music" },
    { "Oldies Music", "Oldies", "Oldies Music" },
    { "Folk Music", "Folk M", "Folk Music" },
    { "Documentary", "Document", "Documentary" },
    { "Alarm Test", "TEST", "Alarm Test" },
    { "Alarm", "Alarm", "Alarm - Alarm !" },
};
