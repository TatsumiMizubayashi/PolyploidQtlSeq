using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// N塩基数の上限値
    /// </summary>
    public class NBaseLimit
    {
        /// <summary>
        /// N塩基数の最小値
        /// </summary>
        public const int MINIMUM = 0;

        /// <summary>
        /// N塩基数の最大値
        /// </summary>
        public const int MAXIMUM = 30;

        /// <summary>
        /// N塩基数の規定値
        /// </summary>
        public const int DEFAULT = 5;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "n";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "nBaseLimit";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Maximum number of N bases. n_base_limit in fastp.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        public const string VALIDATION_ERROR_MESSAGE = "The -n option must be an integer greater than or equal to 0 and less than or equal to 30.";

        /// <summary>
        /// N塩基数の最大値を作成する。
        /// </summary>
        /// <param name="limit">N塩基最大数</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public NBaseLimit(int limit, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, limit, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// N塩基数の最大値を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        internal string ToFastpArg()
        {
            return $"-n {Value}";
        }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }

    }
}
