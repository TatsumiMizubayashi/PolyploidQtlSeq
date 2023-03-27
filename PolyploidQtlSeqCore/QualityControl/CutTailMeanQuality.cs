﻿using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 3'末端トリム時の平均クオリティ
    /// </summary>
    public class CutTailMeanQuality
    {
        /// <summary>
        /// 3'末端トリム時の平均クオリティの最小値
        /// </summary>
        public const int MINIMUM = 1;

        /// <summary>
        /// 3'末端トリム時の平均クオリティの最大値
        /// </summary>
        public const int MAXIMUM = 36;

        /// <summary>
        /// 3'末端トリム時の平均クオリティの規定値
        /// </summary>
        public const int DEFAULT = 20;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "Q";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "cutTailMeanQuality";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Average quality threshold at end trim. cut_tail_mean_quality in fastp.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        public const string VALIDATION_ERROR_MESSAGE = "The -Q option must be an integer greater than 1 and less than or equal to 36.";

        /// <summary>
        /// 3'末端トリム時の平均クオリティを作成する。
        /// </summary>
        /// <param name="quality">平均クオリティ</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public CutTailMeanQuality(int quality, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, quality, parameterDictionary, userOptionDictionary);

            if (quality < MINIMUM || quality > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// 3'末端トリム時の平均クオリティを取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        internal string ToFastpArg()
        {
            return $"--cut_tail_mean_quality {Value}";
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
