using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEffデータベース
    /// </summary>
    public class SnpEffDatabase
    {
        /// <summary>
        /// SnpEffデータベースの規定値
        /// </summary>
        public const string DEFAULT = "";

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "sd";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "snpEffDatabase";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "SnpEff database name.";

        /// <summary>
        /// SnpEffデータベースを作成する。
        /// </summary>
        /// <param name="databaseName">データベース名</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public SnpEffDatabase(string databaseName, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, databaseName, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// SnpEffデータベース名を取得する。
        /// </summary>
        internal string Value { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }
    }
}
