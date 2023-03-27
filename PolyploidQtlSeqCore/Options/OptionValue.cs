namespace PolyploidQtlSeqCore.Options
{
    /// <summary>
    /// オプションの値
    /// </summary>
    internal static class OptionValue
    {
        /// <summary>
        /// オプションの値を変更する必要があるかどうかを調査する。
        /// </summary>
        /// <param name="longName">オプションのLongName</param>
        /// <param name="parameterDictionary">パラメーター辞書</param>
        /// <param name="userSpecifiedLongNameDictionary">ユーザー指定オプションLongNaem辞書</param>
        /// <returns></returns>
        public static bool NeedChangeValue(string longName,
            IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userSpecifiedLongNameDictionary)
        {
            // ユーザー指定されている場合は書き換えの必要なし
            if (userSpecifiedLongNameDictionary.ContainsKey(longName)) return false;

            return parameterDictionary.ContainsKey(longName);
        }

        /// <summary>
        /// 使用する値を取得する。
        /// </summary>
        /// <param name="longName">LongName</param>
        /// <param name="value">値</param>
        /// <param name="parameterDictionary">パラメータファイルの情報</param>
        /// <param name="userOptionDictionary">ユーザー指定オプション情報</param>
        /// <returns></returns>
        public static string GetValue(string longName, string value, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            if (!NeedChangeValue(longName, parameterDictionary, userOptionDictionary)) return value;

            return parameterDictionary[longName];
        }

        /// <summary>
        /// 使用する値を取得する。
        /// </summary>
        /// <param name="longName">longName</param>
        /// <param name="value">値</param>
        /// <param name="parameterDictionary">パラメータファイルの情報</param>
        /// <param name="userOptionDictionary">ユーザー指定オプション情報</param>
        /// <returns></returns>
        public static int GetValue(string longName, int value, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            if (!NeedChangeValue(longName, parameterDictionary, userOptionDictionary)) return value;

            if (!int.TryParse(parameterDictionary[longName], out var newValue))
                throw new InvalidCastException($"Cannot cast {longName} to numeric.");

            return newValue;
        }

        /// <summary>
        /// 使用する値を取得する。
        /// </summary>
        /// <param name="longName">LongName</param>
        /// <param name="value">値</param>
        /// <param name="parameterDictionary">パラメータファイルの情報</param>
        /// <param name="userOptionDictionary">ユーザー指定オプション情報</param>
        /// <returns></returns>
        public static double GetValue(string longName, double value, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            if (!NeedChangeValue(longName, parameterDictionary, userOptionDictionary)) return value;

            if (!double.TryParse(parameterDictionary[longName], out var newValue))
                throw new InvalidCastException($"Cannot cast {longName} to numeric.");

            return newValue;
        }
    }
}
