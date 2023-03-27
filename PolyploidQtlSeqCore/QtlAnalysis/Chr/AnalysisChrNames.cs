using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 解析対象染色体名
    /// </summary>
    public class AnalysisChrNames
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "cn";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "chrNames";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Name of chromosome to be analyzed. If there is more than one, enter them separated by commas. If chrNames is specified, chrSize is disabled.";

        private static readonly char[] _delimiter = new[] { ',' };

        /// <summary>
        /// 解析対象染色体名を作成する。
        /// </summary>
        /// <param name="value">カンマ区切り染色体名</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public AnalysisChrNames(string value, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, value, parameterDictionary, userOptionDictionary);
            Names = string.IsNullOrEmpty(Value)
                ? Array.Empty<string>()
                : Value.Split(_delimiter);
            HasNames = Names.Any();
        }

        /// <summary>
        /// オプション指定された値を取得する。
        /// </summary>
        internal string Value { get; }

        /// <summary>
        /// 解析対象染色体名を取得する。
        /// </summary>
        internal string[] Names { get; }

        /// <summary>
        /// 解析対象染色体名を持っているかどうかを取得する。
        /// </summary>
        internal bool HasNames { get; }

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
