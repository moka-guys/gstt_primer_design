-- Insert default user
INSERT INTO app_users (username, password_hash)
VALUES ('teddy', 'pbkdf2:sha256:260000$GQJRDoO0RGk7kqnM$4fbf462e7949b4064330a0ad3a079e832f320499aa6052e27cec0e2c01d21c5e')
ON CONFLICT (username) DO NOTHING;
