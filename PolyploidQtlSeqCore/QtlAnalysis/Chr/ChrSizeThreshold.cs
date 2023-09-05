using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 染色体長のしきい値
    /// </summary>
    public class ChrSizeThreshold
    {
        /// <summary>
        /// しきい値の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// しきい値の最大値
        /// </summary>
        private const int MAXIMUM = 100_000_000;

        /// <summary>
        /// しきい値の規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 10_000_000;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "cs";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "chrSize";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Threshold for length of chromosomes to be analyzed. Chromosomes with a length more than this value are analyzed.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -cs option must be an integer greater than 1 and less than or equal to 100,000,000.";

        /// <summary>
        /// 染色体長のしきい値
        /// </summary>
        /// <param name="threshold">しきい値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public ChrSizeThreshold(int threshold, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, threshold, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// 染色体長のしきい値
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public ChrSizeThreshold(int threshold)
        {
            if (threshold < MINIMUM || threshold > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(threshold));

            Value = threshold;
        }

        /// <summary>
        /// 染色体サイズのしきい値を取得する。
        /// </summary>
        internal int Value { get; }

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
