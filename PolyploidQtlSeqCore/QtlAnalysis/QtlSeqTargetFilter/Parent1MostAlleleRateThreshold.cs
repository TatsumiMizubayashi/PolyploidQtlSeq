using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// Parent1 MostAlleleRateしきい値
    /// </summary>
    public class Parent1MostAlleleRateThreshold
    {
        /// <summary>
        /// Parent1 MaxAllelRateしきい値の最小値
        /// </summary>
        private const double MINIMUM = 0;

        /// <summary>
        /// Parent1 MaxAllelRateしきい値の最大値
        /// </summary>
        private const double MAXIMUM = 1.0;

        /// <summary>
        /// Parent1 MaxAllelRateしきい値の規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const double DEFAULT = 0.99;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "p1r";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "p1MostAlleleRate";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Most allele frequency for Parent1. Variants exceeding this threshold is considered homozygous.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -p1r option must be a number greater than or equal to 0.0 and less than or equal to 1.0.";

        /// <summary>
        /// Parent1 MostAlleleRateしきい値を作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public Parent1MostAlleleRateThreshold(double threshold, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, threshold, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// Parent1 MostAlleleRateしきい値を作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public Parent1MostAlleleRateThreshold(double threshold)
        {
            if (threshold < MINIMUM || threshold > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(threshold));

            Value = threshold;
        }

        /// <summary>
        /// Parent1 MaxAlleleRateしきい値を取得する。
        /// </summary>
        internal double Value { get; }

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
